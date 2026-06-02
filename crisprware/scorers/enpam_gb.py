"""enPAM+GB Cas12a on-target scorer (Luo et al. 2020).

This is a self-contained port of the inference path used by
`sgrna_modeler.models.SklearnSgrnaModel.predict_seqs` when called with the
enPAM+GB joblib and the `cas12a` enzyme spec. The feature engineering is
vendored verbatim from `sgrna_modeler/features.py` (MIT, Peter DeWeirdt) so
the runtime env needs only scikit-learn + numpy + pandas + biopython.

Run `tools/extract_enpam_gb.py` once to populate
`crisprware/scorers/weights/enpam_gb.joblib` before importing this module.

Architecture: scikit-learn `GradientBoostingRegressor`, RS2 features
(positional-independent 1-/2-mers, positional-dependent 1-/2-mers, GC,
Tm) on the 34-nt Cas12a context with `guide_start=9` (1-indexed),
`guide_length=23`.
"""

from __future__ import annotations

import os
from importlib import resources
from typing import Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

CONTEXT_LEN = 34
GUIDE_START = 9
GUIDE_LENGTH = 23
NTS = ("A", "C", "T", "G")
RS2_FEATURES = (
    "Pos. Ind. 1mer",
    "Pos. Ind. 2mer",
    "Pos. Dep. 1mer",
    "Pos. Dep. 2mer",
    "GC content",
    "Tm",
)
_DEFAULT_WEIGHTS = "enpam_gb.joblib"


def _weights_path() -> str:
    try:
        return str(resources.files(__package__).joinpath("weights", _DEFAULT_WEIGHTS))
    except (AttributeError, ModuleNotFoundError):
        return os.path.join(os.path.dirname(__file__), "weights", _DEFAULT_WEIGHTS)


def is_valid_seq(s: object) -> bool:
    if not isinstance(s, str) or len(s) != CONTEXT_LEN:
        return False
    return all(ch in "ACGTacgt" for ch in s)


def _get_guide_sequence(context: str) -> str:
    return context[(GUIDE_START - 1) : (GUIDE_START - 1 + GUIDE_LENGTH)]


def _gc_fraction(guide: str) -> float:
    return (guide.count("G") + guide.count("C")) / len(guide)


def _one_nt_counts(guide: str) -> dict:
    return {nt: guide.count(nt) / len(guide) for nt in NTS}


def _two_nt_counts(guide: str) -> dict:
    out = {}
    for a in NTS:
        for b in NTS:
            two = a + b
            out[two] = guide.count(two) / (len(guide) - 1)
    return out


def _one_nt_pos(context: str) -> dict:
    out = {}
    for i, ch in enumerate(context):
        pos_key = str(i + 1)
        for nt in NTS:
            out[pos_key + nt] = 1 if ch == nt else 0
    return out


def _two_nt_pos(context: str) -> dict:
    out = {}
    for i in range(len(context) - 1):
        pos_key = str(i + 1)
        pair = context[i : i + 2]
        for a in NTS:
            for b in NTS:
                key = pos_key + a + b
                out[key] = 1 if pair == a + b else 0
    return out


def _thermo(guide: str, context: str) -> dict:
    from Bio.SeqUtils import MeltingTemp

    third = len(guide) // 3
    return {
        "Tm, context": MeltingTemp.Tm_NN(context),
        "Tm, start": MeltingTemp.Tm_NN(guide[:third]),
        "Tm, mid": MeltingTemp.Tm_NN(guide[third : 2 * third]),
        "Tm, end": MeltingTemp.Tm_NN(guide[2 * third :]),
    }


# Vectorized featurizer helpers (much faster than per-sequence Python loops).
# Lookup table: byte value (ASCII) -> nt index in NTS order (A=0, C=1, T=2, G=3).
_NT_LUT = np.zeros(256, dtype=np.int8)
for _c, _i in zip(b"ACTG", range(4)):
    _NT_LUT[_c] = _i
    _NT_LUT[_c | 0x20] = _i  # also accept lowercase
_EYE4 = np.eye(4, dtype=np.float64)
_EYE16 = np.eye(16, dtype=np.float64)


def _featurize_array(seqs: Sequence[str]) -> np.ndarray:
    """Vectorized featurizer returning a (N, 689) float64 array.

    Column order matches the original `featurize()` exactly so the model's
    learned column ordering still applies. Per-sequence Python loops are
    replaced with numpy ops for the 685 positional / count / GC features;
    the 4 Tm features still need a Python loop (Biopython's Tm_NN is
    per-sequence) but those are <40% of total cost.
    """
    n = len(seqs)
    # Encode N sequences as (N, 34) int8 of nt indices 0-3
    joined = "".join(s.upper() for s in seqs)
    raw = np.frombuffer(joined.encode("ascii"), dtype=np.uint8).reshape(n, CONTEXT_LEN)
    arr = _NT_LUT[raw].astype(np.int64)  # (N, 34)
    guide = arr[:, GUIDE_START - 1 : GUIDE_START - 1 + GUIDE_LENGTH]  # (N, 23)

    # 1-mer counts on the 23-nt guide; one-hot then sum along position.
    counts_1mer = _EYE4[guide].sum(axis=1)  # (N, 4)  in NTS order ACTG
    frac_1mer = counts_1mer / GUIDE_LENGTH
    # GC content = (C + G) / 23  (NTS indices: C=1, G=3)
    gc = (counts_1mer[:, 1] + counts_1mer[:, 3]) / GUIDE_LENGTH

    # 2-mer counts on the guide. Upstream uses Python's str.count, which counts
    # *non-overlapping* matches (e.g. "AAAA".count("AA") == 2, not 3). The naive
    # positional-overlap count diverges for homo-dinucs (AA, CC, TT, GG) where
    # a run of length k yields floor(k/2) non-overlapping matches, not k-1.
    # 16 str.count() calls per sequence is cheap (~16μs/seq); leave it scalar.
    frac_2mer = np.empty((n, 16), dtype=np.float64)
    dinucs = [a + b for a in NTS for b in NTS]
    denom = float(GUIDE_LENGTH - 1)
    for i, s in enumerate(seqs):
        g = s.upper()[GUIDE_START - 1 : GUIDE_START - 1 + GUIDE_LENGTH]
        for j, dn in enumerate(dinucs):
            frac_2mer[i, j] = g.count(dn) / denom

    # Positional 1-mer over the full 34-nt context: one-hot reshape (position outer, nt inner).
    pos_1mer = _EYE4[arr].reshape(n, CONTEXT_LEN * 4)  # (N, 136)

    # Positional 2-mer over the 34-nt context: 33 positions × 16 dinucs.
    ctx_dinuc = arr[:, :-1] * 4 + arr[:, 1:]  # (N, 33)
    pos_2mer = _EYE16[ctx_dinuc].reshape(n, (CONTEXT_LEN - 1) * 16)  # (N, 528)

    # Tm features still need per-sequence Biopython calls.
    from Bio.SeqUtils import MeltingTemp as MT

    tm = np.empty((n, 4), dtype=np.float64)
    for i, s in enumerate(seqs):
        s_up = s.upper()
        g = s_up[GUIDE_START - 1 : GUIDE_START - 1 + GUIDE_LENGTH]
        third = GUIDE_LENGTH // 3
        tm[i, 0] = float(MT.Tm_NN(s_up))
        tm[i, 1] = float(MT.Tm_NN(g[:third]))
        tm[i, 2] = float(MT.Tm_NN(g[third : 2 * third]))
        tm[i, 3] = float(MT.Tm_NN(g[2 * third :]))

    # Concatenate in the original column order
    return np.concatenate(
        [
            gc.reshape(-1, 1),  # 1
            frac_1mer,  # 4
            frac_2mer,  # 16
            pos_1mer,  # 136
            pos_2mer,  # 528
            tm,  # 4
        ],
        axis=1,
    )  # (N, 689)


def _feature_names() -> list:
    """Generate the 689 feature names in the same order as the legacy dict."""
    names = ["GC content"]
    names += list(NTS)  # "A","C","T","G"
    names += [a + b for a in NTS for b in NTS]  # AA, AC, ..., GG  (16)
    names += [f"{i + 1}{nt}" for i in range(CONTEXT_LEN) for nt in NTS]  # 1A, 1C, ..., 34G
    names += [f"{i + 1}{a}{b}" for i in range(CONTEXT_LEN - 1) for a in NTS for b in NTS]
    names += ["Tm, context", "Tm, start", "Tm, mid", "Tm, end"]
    return names


def featurize(seqs: Sequence[str]) -> pd.DataFrame:
    """RS2-style features for Cas12a 34-nt context.

    Wraps `_featurize_array` so the public API still returns a labeled
    DataFrame; the hot path in `predict()` calls `_featurize_array`
    directly to skip pandas overhead.
    """
    return pd.DataFrame(_featurize_array(seqs), columns=_feature_names())


def _install_sklearn_legacy_aliases() -> None:
    """Map sklearn 0.21.x module paths onto current sklearn so the upstream
    enPAM+GB pickle (pickled with sklearn 0.21.2) can unpickle in sklearn 1.x.

    Predictions are bit-exact (mean abs diff < 1e-8) once these aliases land —
    the tree data structures didn't change, only the module layout did.
    """
    import sys

    try:
        import sklearn.ensemble._gb

        sys.modules.setdefault("sklearn.ensemble.gradient_boosting", sklearn.ensemble._gb)
    except ImportError:
        pass
    try:
        import sklearn.tree._classes

        sys.modules.setdefault("sklearn.tree.tree", sklearn.tree._classes)
    except ImportError:
        pass


def load_model(weights_path: Optional[str] = None):
    """Load the GradientBoostingRegressor from the vendored joblib."""
    import warnings
    import joblib

    path = weights_path or _weights_path()
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"enPAM+GB weights not found at {path}. Run `python tools/extract_enpam_gb.py` once to populate it."
        )

    _install_sklearn_legacy_aliases()
    with warnings.catch_warnings():
        # Suppress the InconsistentVersionWarning that sklearn 1.x raises when
        # loading a 0.21 pickle. Bit-exact equivalence is gated by tests.
        warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")
        try:
            from sklearn.exceptions import InconsistentVersionWarning

            warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
        except ImportError:
            pass
        return joblib.load(path)


def predict(
    seqs: Sequence[str],
    model=None,
    batch_size: int = 4096,
    weights_path: Optional[str] = None,
) -> List[float]:
    """Score sequences. Invalid (non-34-nt-ACGT) sequences return NaN."""
    if model is None:
        model = load_model(weights_path)

    n = len(seqs)
    out = [float("nan")] * n
    valid_idx = [i for i, s in enumerate(seqs) if is_valid_seq(s)]
    if not valid_idx:
        return out

    valid_seqs = [seqs[i].upper() for i in valid_idx]
    for start in range(0, len(valid_seqs), batch_size):
        chunk = valid_seqs[start : start + batch_size]
        # Skip the DataFrame wrapping entirely in the hot path -- direct
        # numpy array also dodges sklearn's "fitted without feature names"
        # warning.
        features = _featurize_array(chunk)
        scores = np.asarray(model.predict(features), dtype=float)
        for slot_local, score in enumerate(scores):
            out[valid_idx[start + slot_local]] = float(score)
    return out


def compute_enpam_gb_scores(
    sequences: Iterable[str],
    threads: int = 1,
    chunk_size: int = 100_000,
    weights_path: Optional[str] = None,
) -> List[float]:
    """Chunked entry mirroring `compute_rs3_scores` in score_guides.py."""
    seqs = list(sequences)
    model = load_model(weights_path)
    scores: List[float] = []
    for i in range(0, len(seqs), chunk_size):
        scores.extend(predict(seqs[i : i + chunk_size], model=model))
    return scores


def score_file(
    input_path: str,
    output_path: str,
    seq_col: str = "context",
    out_col: str = "enpam_gb_score",
    batch_size: int = 4096,
) -> None:
    """Read a TSV, score `seq_col`, append `out_col`, write TSV."""
    df = pd.read_csv(input_path, sep="\t", dtype=str)
    if df.columns[0].startswith("#"):
        df.rename(columns={df.columns[0]: df.columns[0].lstrip("#")}, inplace=True)
    if seq_col not in df.columns:
        raise ValueError(f"Column {seq_col!r} not found. Available: {list(df.columns)}")

    seqs = df[seq_col].fillna("").tolist()
    scores = predict(seqs, batch_size=batch_size)
    df[out_col] = [None if (s != s) else round(float(s), 8) for s in scores]
    df.to_csv(output_path, sep="\t", index=False)
