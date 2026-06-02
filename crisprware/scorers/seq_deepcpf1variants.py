"""seq-DeepCpf1variants: per-variant Cas12a on-target scorers (Chen et al. 2025).

23 trained models -- one per Cas12a variant in Chen 2025's library, all
sharing the EnCas12a sequence-only architecture (identical to enseq-DeepCpf1
in `crisprware.scorers.enseq_deepcpf1`). The architecture and weight-loading
logic is reused from that module; this module's only job is to map a variant
name to the right weights file.

Notable: this is the only public scorer for engineered variants like
`AsCas12a Ultra` (= 2xNLS-Cas12a), `enAsCas12a-HF1`, `HyperFi-AsCas12a`,
and the Lb / Fn / Ce / Eb orthologs. Per Chen Fig 4b, each variant model
gets Spearman 0.79-0.91 on its own held-out test data.

Variants (matching the .npz filenames in weights/chen_2025/variants/):

    AsCas12a, AsCas12a-Plus, AsCas12aRR, AsCas12aRVR, AsCas12a_Ultra,
    enAsCas12a-HF1, HyperFi-AsCas12a,
    LbCas12a, LbCas12a-Plus, LbCas12aK538R, LbCas12aRR, LbCas12aRVR, LbCas12aRVRR,
    Lb2Cas12a, Lb2Cas12aK518R,
    FnCas12a, FnCas12aRVR, eaFnCas12a,
    CeCas12a, EbCas12a, enEbCas12a,
    iCas12a_mut2C-W, iCas12a_mut2C-WF.
"""

from __future__ import annotations

import os
from importlib import resources
from typing import Iterable, List, Optional, Sequence


from crisprware.scorers import enseq_deepcpf1
from crisprware.scorers.enseq_deepcpf1 import (
    CONTEXT_LEN,
)

VARIANTS = (
    "AsCas12a",
    "AsCas12a-Plus",
    "AsCas12aRR",
    "AsCas12aRVR",
    "AsCas12a_Ultra",
    "enAsCas12a-HF1",
    "HyperFi-AsCas12a",
    "LbCas12a",
    "LbCas12a-Plus",
    "LbCas12aK538R",
    "LbCas12aRR",
    "LbCas12aRVR",
    "LbCas12aRVRR",
    "Lb2Cas12a",
    "Lb2Cas12aK518R",
    "FnCas12a",
    "FnCas12aRVR",
    "eaFnCas12a",
    "CeCas12a",
    "EbCas12a",
    "enEbCas12a",
    "iCas12a_mut2C-W",
    "iCas12a_mut2C-WF",
)


def _normalize_variant(name: str) -> str:
    """Map user-supplied variant name to the canonical filename stem.

    Accepts case-insensitive input and allows ' '/'-'/'_' / '('/')' interchange.
    """
    norm = name.strip().replace(" ", "_").replace("(", "_").replace(")", "")
    canon_map = {v.lower(): v for v in VARIANTS}
    if norm.lower() in canon_map:
        return canon_map[norm.lower()]
    # Also try without dashes (e.g., "ascas12aplus" -> "AsCas12a-Plus")
    no_dash = norm.replace("-", "").lower()
    for v in VARIANTS:
        if v.replace("-", "").lower() == no_dash:
            return v
    raise ValueError(f"Unknown Cas12a variant {name!r}. Available: {', '.join(VARIANTS)}")


def _weights_path(variant: str) -> str:
    canonical = _normalize_variant(variant)
    sub = os.path.join("chen_2025", "variants", f"{canonical}.npz")
    try:
        return str(resources.files(enseq_deepcpf1.__package__).joinpath("weights", sub))
    except (AttributeError, ModuleNotFoundError):
        return os.path.join(os.path.dirname(enseq_deepcpf1.__file__), "weights", sub)


def load_model(variant: str, weights_path: Optional[str] = None):
    """Build the EnCas12a architecture and load the variant-specific weights."""
    path = weights_path or _weights_path(variant)
    return enseq_deepcpf1.load_model(weights_path=path)


def predict(
    seqs: Sequence[str],
    variant: str,
    model=None,
    batch_size: int = 1024,
    weights_path: Optional[str] = None,
    input_length: int = CONTEXT_LEN,
) -> List[float]:
    """Score sequences with the variant-specific seq-DeepCpf1variants model."""
    if model is None:
        model = load_model(variant, weights_path=weights_path)
    return enseq_deepcpf1.predict(seqs, model=model, batch_size=batch_size, input_length=input_length)


def compute_variant_scores(
    sequences: Iterable[str],
    variant: str,
    threads: int = 1,
    chunk_size: int = 100_000,
    weights_path: Optional[str] = None,
) -> List[float]:
    """Chunked entry mirroring `compute_enseq_deepcpf1_scores`."""
    seqs = list(sequences)
    model = load_model(variant, weights_path=weights_path)
    scores: List[float] = []
    for i in range(0, len(seqs), chunk_size):
        scores.extend(predict(seqs[i : i + chunk_size], variant=variant, model=model))
    return scores


def score_file(
    input_path: str,
    output_path: str,
    variant: str,
    seq_col: str = "context",
    out_col: Optional[str] = None,
    batch_size: int = 1024,
) -> None:
    """Read a TSV, score `seq_col` with the variant model, append `out_col`."""
    import pandas as pd

    canonical = _normalize_variant(variant)
    if out_col is None:
        out_col = f"seq_deepcpf1variants_{canonical.replace('-', '_')}_score"

    df = pd.read_csv(input_path, sep="\t", dtype=str)
    if df.columns[0].startswith("#"):
        df.rename(columns={df.columns[0]: df.columns[0].lstrip("#")}, inplace=True)
    if seq_col not in df.columns:
        raise ValueError(f"Column {seq_col!r} not found.")
    seqs = df[seq_col].fillna("").tolist()
    scores = predict(seqs, variant=canonical, batch_size=batch_size)
    df[out_col] = [None if (s != s) else round(float(s), 8) for s in scores]
    df.to_csv(output_path, sep="\t", index=False)
