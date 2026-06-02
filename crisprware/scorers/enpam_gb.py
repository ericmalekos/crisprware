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


def featurize(seqs: Sequence[str]) -> pd.DataFrame:
    """RS2-style features for Cas12a 34-nt context.

    Schema (column order) mirrors `sgrna_modeler.features.featurize_guides`
    with the default RS2 feature list. The order matters because the
    upstream model was trained with these features iterated through Python
    dict insertion order; we replicate that exactly.
    """
    rows: List[dict] = []
    for ctx in seqs:
        guide = _get_guide_sequence(ctx)
        d: dict = {}
        d.update({"GC content": _gc_fraction(guide)})
        d.update(_one_nt_counts(guide))
        d.update(_two_nt_counts(guide))
        d.update(_one_nt_pos(ctx))
        d.update(_two_nt_pos(ctx))
        d.update(_thermo(guide, ctx))
        rows.append(d)
    return pd.DataFrame(rows)


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
        features = featurize(chunk)
        # Use .to_numpy() to silence sklearn's "fitted without feature names"
        # warning — the upstream model was trained on a numpy array, not a DataFrame.
        scores = np.asarray(model.predict(features.to_numpy()), dtype=float)
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
