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


def _build_tm_tables():
    """Build per-dinucleotide enthalpy/entropy tables from Biopython's
    DNA_NN3 (Allawi & SantaLucia 1997). Returns two length-16 numpy arrays
    indexed by dinuc_idx = nt1*4 + nt2 in NTS order (A=0, C=1, T=2, G=3).

    DNA_NN3 has 10 unique dinucleotide entries; the other 6 are obtained
    by reverse-complement symmetry (mirrors Biopython's neighbors-reversed
    fallback in `Tm_NN`'s zip loop).
    """
    from Bio.SeqUtils.MeltingTemp import DNA_NN3

    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    h = np.zeros(16, dtype=np.float64)
    s = np.zeros(16, dtype=np.float64)
    for i1, n1 in enumerate(NTS):  # NTS = ("A","C","T","G")
        for i2, n2 in enumerate(NTS):
            idx = i1 * 4 + i2
            key = f"{n1}{n2}/{comp[n1]}{comp[n2]}"
            if key in DNA_NN3:
                h[idx], s[idx] = DNA_NN3[key]
            elif key[::-1] in DNA_NN3:
                h[idx], s[idx] = DNA_NN3[key[::-1]]
            else:
                raise ValueError(f"missing NN entry for {key}")
    return h, s


_TM_H_TABLE, _TM_S_TABLE = _build_tm_tables()


def _tm_nn_fast(arr_int: np.ndarray, dnac1: float = 25.0, dnac2: float = 25.0, Na: float = 50.0) -> np.ndarray:
    """Vectorized Biopython Tm_NN equivalent for batched fixed-length DNA.

    Matches `MeltingTemp.Tm_NN(seq, dnac1=25, dnac2=25, Na=50)` to within
    float64 noise on canonical ACTG inputs. Skips the per-call input
    validation / `Bio.Seq` wrapping / mismatch-table-lookup overhead that
    accounts for >60% of Biopython's per-sequence cost.

    Args:
        arr_int: (N, L) int array, values in {0,1,2,3} encoding A/C/T/G
            in NTS order.
    Returns:
        (N,) float64 array of Tm values in °C.
    """
    n, L = arr_int.shape
    # DNA_NN3 "init", "init_oneG/C", "init_allA/T", "init_5T/A" all (0,0),
    # so only init_A/T (2.3, 4.1) and init_G/C (0.1, -2.8) terminal contributions matter.
    # NTS=("A","C","T","G") -> A/T are indices 0, 2; G/C are 1, 3.
    at_first = ((arr_int[:, 0] == 0) | (arr_int[:, 0] == 2)).astype(np.float64)
    at_last = ((arr_int[:, -1] == 0) | (arr_int[:, -1] == 2)).astype(np.float64)
    at_count = at_first + at_last
    gc_count = 2.0 - at_count
    delta_h = 2.3 * at_count + 0.1 * gc_count
    delta_s = 4.1 * at_count + (-2.8) * gc_count

    # Zip: sum dinucleotide H/S contributions.
    dinuc = arr_int[:, :-1] * 4 + arr_int[:, 1:]
    delta_h = delta_h + _TM_H_TABLE[dinuc].sum(axis=1)
    delta_s = delta_s + _TM_S_TABLE[dinuc].sum(axis=1)

    # Salt correction (saltcorr=5 default in Tm_NN): 0.368 * (L-1) * ln([Mon])
    # added to delta_s. With Na=50 only, [Mon] = 50e-3 M.
    mon = Na * 1e-3
    delta_s = delta_s + 0.368 * (L - 1) * np.log(mon)

    # Tm = (1000 * dH) / (dS + R * ln(k)) - 273.15
    #   k = (dnac1 - dnac2/2) * 1e-9
    k = (dnac1 - dnac2 / 2.0) * 1e-9
    R = 1.987
    return (1000.0 * delta_h) / (delta_s + R * np.log(k)) - 273.15


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

    # Tm features via the vectorized Tm_NN equivalent. Four batched calls
    # replace 4N per-sequence Biopython invocations (~6-10x faster overall).
    third = GUIDE_LENGTH // 3  # = 7 for 23-nt guide
    tm = np.empty((n, 4), dtype=np.float64)
    tm[:, 0] = _tm_nn_fast(arr)  # full 34-nt context
    tm[:, 1] = _tm_nn_fast(guide[:, :third])  # 7-mer start
    tm[:, 2] = _tm_nn_fast(guide[:, third : 2 * third])  # 7-mer middle
    tm[:, 3] = _tm_nn_fast(guide[:, 2 * third :])  # 9-mer end

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
