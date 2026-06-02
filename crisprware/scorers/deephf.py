"""DeepHF on-target activity scorer (Wang et al. 2019, Nat Commun).

TF2/Keras port of the BiLSTM + bio-features model from Wang, Zhang, Wang
et al. 2019 (DOI 10.1038/s41467-019-12281-8). Three variants (matching
the three .hdf5 files preserved via Bioconductor's crisprScoreData /
ExperimentHub) are supported:

    wt_u6  -- wildtype SpCas9, U6 promoter        (DeepWt_U6.hdf5)
    esp    -- enhanced SpCas9 (eSpCas9-1.1)       (esp_rnn_model.hdf5)
    hf     -- high-fidelity SpCas9 (SpCas9-HF1)   (hf_rnn_model.hdf5)

The fourth published variant (`wt_t7`, T7 promoter) is not present in the
crisprScoreData cache that ships with this repo; if needed, users can
download it from ExperimentHub resource EH6128 or contact the authors.

Original GitHub repo (`izhangcd/DeepHF`) has been deleted; weights are
preserved by Bioconductor crisprScoreData. The Wang 2019 article and
its data are open-access on Nature Comm.

Input: 23-nt protospacer + PAM (20 spacer + 3 PAM, NGG). The first 21
characters (20 spacer + 1 PAM nt) are integer-tokenized to a 22-long
vector (START token + 21 chars, padded with zeros at the front by
Keras's pad_sequences; in practice the first slot is the START=1 token
and the rest are the encoded characters A=2, T=3, C=4, G=5).

Bio features (11 dims, see feature_util.py upstream):
    stem, dG, dG_binding_20, dg_binding_7to20,  # 4 RNA folding features
    GC>10, GC<10, GC count,                     # 3 GC features
    Tm global, 5mer end, 8mer middle, 4mer start  # 4 melting temp features

ViennaRNA (Python `RNA` module, pip-installable) is required only for
the `dG` and `stem` features; the other 9 features are pure
numpy/Biopython.
"""

from __future__ import annotations

import os
from importlib import resources
from typing import Iterable, List, Optional, Sequence

import numpy as np

# Token IDs match prediction_util.py's make_data: START=1, A=2, T=3, C=4, G=5.
SEQ_LEN = 23
TOKEN_LEN = 22
_CHAR_TO_ID = {"A": 2, "T": 3, "C": 4, "G": 5, "a": 2, "t": 3, "c": 4, "g": 5}
_START_ID = 1

# scaffold sequence appended to the 20-bp guide for RNAfold (from feature_util.py:336)
_TRACR_SCAFFOLD = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTT"

# Nearest-neighbor dinucleotide free energies for dG_binding (feature_util.py:418-422)
_DG_DINUC = {
    "aa": -0.2,
    "tt": -1.0,
    "at": -0.9,
    "ta": -0.6,
    "ca": -1.6,
    "tg": -0.9,
    "ct": -1.8,
    "ag": -0.9,
    "ga": -1.5,
    "tc": -1.3,
    "gt": -2.1,
    "ac": -1.1,
    "cg": -1.7,
    "gc": -2.7,
    "gg": -2.1,
    "cc": -2.9,
}
_DG_INIT = 3.1

# Expected stem pattern at positions [18:50] of the 99-bp folded structure
_STEM_PATTERN = "(((((((((.((((....))))...)))))))"

# Tm segment slices on the 21-mer (Tm_feature, feature_util.py:150)
_TM_SEGMENTS = [(15, 21), (4, 13), (0, 4)]

_VARIANT_TO_WEIGHTS = {
    "wt_u6": "DeepWt_U6.hdf5",
    "esp": "esp_rnn_model.hdf5",
    "hf": "hf_rnn_model.hdf5",
}


def _weights_path(variant: str) -> str:
    fn = _VARIANT_TO_WEIGHTS.get(variant)
    if fn is None:
        raise ValueError(f"Unknown DeepHF variant {variant!r}. Available: {sorted(_VARIANT_TO_WEIGHTS)}")
    try:
        return str(resources.files(__package__).joinpath("weights", "wang_2019_deephf", fn))
    except (AttributeError, ModuleNotFoundError):
        return os.path.join(os.path.dirname(__file__), "weights", "wang_2019_deephf", fn)


def is_valid_seq(s: object, length: int = SEQ_LEN) -> bool:
    if not isinstance(s, str) or len(s) != length:
        return False
    return all(ch in _CHAR_TO_ID for ch in s)


def encode_sequence(seqs_23: Sequence[str]) -> np.ndarray:
    """Integer-encode 23-nt sequences for DeepHF's BiLSTM input.

    Replicates `make_data` from upstream prediction_util.py:
      - take `seq[:21]` (20 spacer + 1 PAM nt)
      - tokenize with {A:2, T:3, C:4, G:5}
      - prepend START=1
      - left-pad to length 22 (no-op since input is exactly 22 chars)
    Output dtype is int32 to match Keras Embedding expectations.
    """
    n = len(seqs_23)
    arr = np.zeros((n, TOKEN_LEN), dtype=np.int32)
    for i, s in enumerate(seqs_23):
        if not isinstance(s, str) or len(s) < 21:
            continue
        arr[i, 0] = _START_ID
        for j, ch in enumerate(s[:21]):
            tok = _CHAR_TO_ID.get(ch)
            if tok is not None:
                arr[i, j + 1] = tok
    return arr


def _dG_binding(seq: str) -> float:
    """Pure-Python dinucleotide free-energy estimator (feature_util.py:416)."""
    s = seq.lower().replace("u", "t")
    total = _DG_INIT
    for i in range(len(s) - 1):
        total += _DG_DINUC.get(s[i : i + 2], 0.0)
    return total


def _gc_features(s21: str) -> tuple:
    """(gc_above_10, gc_below_10, gc_count) computed over s21[:20] (feature_util.py:134, 171)."""
    g = s21[:20].count("G") + s21[:20].count("C")
    return (1 if g > 10 else 0, 1 if g < 10 else 0, g)


def _tm_features(s21: str) -> tuple:
    """4 melting-temperature features (feature_util.py:147).

    Upstream uses `Bio.SeqUtils.MeltingTemp.Tm_staluc(seq, rna=False)`, which
    Biopython deprecated in 1.78+ and removed entirely later. The bundled
    Tm_NN with these exact args is the documented equivalent: per Biopython's
    own docstring, `Tm_staluc(s, dnac=50, saltc=50, rna=0)` is just
    `Tm_NN(s, dnac1=25, dnac2=25, Na=50)`.
    """
    from Bio.SeqUtils import MeltingTemp as Tm

    def _tm(seq):
        return float(Tm.Tm_NN(seq, dnac1=25, dnac2=25, Na=50))

    seg_full = _tm(s21)
    seg_5mer_end = _tm(s21[_TM_SEGMENTS[0][0] : _TM_SEGMENTS[0][1]])
    seg_8mer_mid = _tm(s21[_TM_SEGMENTS[1][0] : _TM_SEGMENTS[1][1]])
    seg_4mer_start = _tm(s21[_TM_SEGMENTS[2][0] : _TM_SEGMENTS[2][1]])
    return (seg_full, seg_5mer_end, seg_8mer_mid, seg_4mer_start)


def _rna_fold_features(s21: str) -> tuple:
    """(stem, dG) via ViennaRNA's RNAfold.

    Upstream (feature_util.py:336-374) folds the 20-bp guide alone AND
    the 20-bp guide + 79-bp tracr scaffold (99 bp total). The free
    energy `dG` comes from the **20-bp fold alone**; the `stem` flag
    comes from comparing the 99-bp fold's structure positions [18:50]
    against an expected stem-loop dot-bracket pattern.
    """
    import RNA

    guide_rna = s21[:20].replace("T", "U")
    fc20 = RNA.fold_compound(guide_rna)
    _, dG = fc20.mfe()

    fc99 = RNA.fold_compound(guide_rna + _TRACR_SCAFFOLD.replace("T", "U"))
    structure99, _ = fc99.mfe()
    aligned_stem = structure99[18 : 18 + len(_STEM_PATTERN)]
    stem = 1 if aligned_stem == _STEM_PATTERN else 0

    return (stem, float(dG))


def compute_bio_features(seqs_23: Sequence[str]) -> np.ndarray:
    """Compute the 11-dim bio feature vector for each 23-nt sequence.

    Column order (must match the upstream model's `bio_input` shape (11,)):
      0: stem            (ViennaRNA)
      1: dG              (ViennaRNA)
      2: dG_binding_20   (pure Python dinuc lookup)
      3: dg_binding_7to20 (pure Python dinuc lookup)
      4: GC > 10         (boolean from countGC)
      5: GC < 10         (boolean from countGC)
      6: GC count        (countGC on s[:20])
      7: Tm global       (Biopython Tm_staluc on full 21mer)
      8: 5mer_end        (Tm_staluc on s[15:21])
      9: 8mer_middle     (Tm_staluc on s[4:13])
     10: 4mer_start      (Tm_staluc on s[0:4])
    """
    n = len(seqs_23)
    feat = np.zeros((n, 11), dtype=np.float32)
    for i, s in enumerate(seqs_23):
        if not isinstance(s, str) or len(s) < 21:
            continue
        s21 = s[:21]
        stem, dG = _rna_fold_features(s21)
        dGb20 = _dG_binding(s21[:20])
        dGb7to20 = _dG_binding(s21[7:20])
        gc_above, gc_below, gc_cnt = _gc_features(s21)
        tm_full, tm_5e, tm_8m, tm_4s = _tm_features(s21)
        feat[i] = [
            stem,
            dG,
            dGb20,
            dGb7to20,
            gc_above,
            gc_below,
            gc_cnt,
            tm_full,
            tm_5e,
            tm_8m,
            tm_4s,
        ]
    return feat


def load_model(variant: str = "wt_u6", weights_path: Optional[str] = None):
    """Load the bundled Keras .hdf5 model for the given variant."""
    import warnings
    import tensorflow as tf

    path = weights_path or _weights_path(variant)
    if not os.path.exists(path):
        raise FileNotFoundError(f"DeepHF weights for variant {variant!r} not found at {path}.")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return tf.keras.models.load_model(path, compile=False)


def predict(
    seqs: Sequence[str],
    variant: str = "wt_u6",
    model=None,
    batch_size: int = 4096,
    weights_path: Optional[str] = None,
) -> List[float]:
    """Score 23-nt protospacer+PAM sequences. Invalid sequences return NaN."""
    if model is None:
        model = load_model(variant, weights_path=weights_path)

    n = len(seqs)
    out = [float("nan")] * n
    valid_idx = [i for i, s in enumerate(seqs) if is_valid_seq(s)]
    if not valid_idx:
        return out

    valid_seqs = [seqs[i].upper() for i in valid_idx]
    x_seq = encode_sequence(valid_seqs)
    x_bio = compute_bio_features(valid_seqs)
    y = model.predict([x_seq, x_bio], batch_size=batch_size, verbose=0).ravel()
    # Upstream clips to [0, 1] (prediction_util.py:103)
    y = np.clip(y, 0.0, 1.0)
    for slot, v in zip(valid_idx, y):
        out[slot] = float(v)
    return out


def compute_deephf_scores(
    sequences: Iterable[str],
    variant: str = "wt_u6",
    threads: int = 1,
    chunk_size: int = 100_000,
    weights_path: Optional[str] = None,
) -> List[float]:
    """Chunked entry mirroring `compute_rs3_scores` in score_guides.py."""
    seqs = list(sequences)
    model = load_model(variant, weights_path=weights_path)
    scores: List[float] = []
    for i in range(0, len(seqs), chunk_size):
        scores.extend(predict(seqs[i : i + chunk_size], variant=variant, model=model))
    return scores


def score_file(
    input_path: str,
    output_path: str,
    variant: str = "wt_u6",
    seq_col: str = "sequence",
    out_col: Optional[str] = None,
    batch_size: int = 4096,
) -> None:
    """Read a TSV, score `seq_col` with DeepHF, append `out_col`."""
    import pandas as pd

    if out_col is None:
        out_col = f"deephf_{variant}_score"
    df = pd.read_csv(input_path, sep="\t", dtype=str)
    if df.columns[0].startswith("#"):
        df.rename(columns={df.columns[0]: df.columns[0].lstrip("#")}, inplace=True)
    if seq_col not in df.columns:
        raise ValueError(f"Column {seq_col!r} not found.")
    seqs = df[seq_col].fillna("").tolist()
    scores = predict(seqs, variant=variant, batch_size=batch_size)
    df[out_col] = [None if (s != s) else round(float(s), 8) for s in scores]
    df.to_csv(output_path, sep="\t", index=False)
