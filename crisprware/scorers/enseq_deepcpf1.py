"""enseq-DeepCpf1 Cas12a on-target scorer (Chen et al. 2025).

TF2 port of the sequence-only CNN from Chen, Wu, Wang, Liu et al. 2025,
Nat Commun 16:3022 (DOI 10.1038/s41467-025-57150-9). Architecture and
trained weights are released under MIT (capsule code) and CC0-1.0 (capsule
data) via Code Ocean capsule 9398276 -- this module ports the PyTorch
EnCas12a model and inference path to TF2/Keras and loads the weights from
a framework-agnostic .npz extracted by `tools/extract_chen_2025_weights.py`.

Architecture
------------
Input:    one-hot (4, 31) channels-first, encoding A=0, T=1, C=2, G=3.
Conv:     Conv1D(128, k=5, padding="same") x 2, ReLU, no pooling.
FC:       Flatten -> Dense(128, ReLU) x 2 -> Dense(1, sigmoid).
Output:   sigmoid in [0, 1].

The 31-nt input layout in Chen's training: 3 nt 5' flanking + 4 nt PAM
(TTTV by default) + 20 nt of the 23-nt AsCas12a protospacer + 4 nt that
overlaps the protospacer tail / 3' flank. crisprware's standard Cas12a
context is 34 nt (4 + 4 + 23 + 3); the right slice is `context[1:32]`,
verified to equal Chen's `seq` column in the capsule's published datasets.
"""

from __future__ import annotations

import os
from importlib import resources
from typing import Iterable, List, Optional, Sequence

import numpy as np

SEQ_LEN = 31
CONTEXT_LEN = 34
NT_INDEX = {
    "A": 0,
    "T": 1,
    "C": 2,
    "G": 3,
    "a": 0,
    "t": 1,
    "c": 2,
    "g": 3,
}

_DEFAULT_WEIGHTS = os.path.join("chen_2025", "enseq_deepcpf1.npz")

_LAYER_KEYS = [
    ("conv1", "conv1_w", "conv1_b"),
    ("conv2", "conv2_w", "conv2_b"),
    ("fc1", "fc1_w", "fc1_b"),
    ("fc2", "fc2_w", "fc2_b"),
    ("out", "out_w", "out_b"),
]


def _weights_path(name: str = _DEFAULT_WEIGHTS) -> str:
    try:
        return str(resources.files(__package__).joinpath("weights", name))
    except (AttributeError, ModuleNotFoundError):
        return os.path.join(os.path.dirname(__file__), "weights", name)


def _build_model():
    """Build the TF2 architecture mirroring EnCas12a from the capsule's model.py.

    Internally uses channels_last (CPU TF doesn't support NCHW Conv) with a
    Permute((2, 1)) before Flatten so the flat output is in channel-major
    order matching PyTorch's `view(B, -1)`. With that ordering the FC1
    weights load directly from the PyTorch state_dict (just transposed) --
    no row permutation needed.
    """
    import tensorflow as tf

    L = tf.keras.layers

    inp = L.Input(shape=(SEQ_LEN, 4), name="input")
    x = L.Conv1D(128, 5, padding="same", activation="relu", name="conv1")(inp)
    x = L.Conv1D(128, 5, padding="same", activation="relu", name="conv2")(x)
    x = L.Permute((2, 1), name="to_channels_first")(x)
    x = L.Flatten(name="flatten")(x)
    x = L.Dense(128, activation="relu", name="fc1")(x)
    x = L.Dense(128, activation="relu", name="fc2")(x)
    out = L.Dense(1, activation="sigmoid", name="out")(x)
    return tf.keras.Model(inp, out, name="enseq_deepcpf1")


def load_model(weights_path: Optional[str] = None):
    """Build the model and load weights from an .npz produced by the extraction tool."""
    model = _build_model()
    path = weights_path or _weights_path()
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"enseq-DeepCpf1 weights not found at {path}. "
            f"Run `python tools/extract_chen_2025_weights.py` once to populate it."
        )
    w = np.load(path)
    for layer_name, w_key, b_key in _LAYER_KEYS:
        model.get_layer(layer_name).set_weights([w[w_key], w[b_key]])
    return model


def is_valid_seq(s: object, length: int = CONTEXT_LEN) -> bool:
    if not isinstance(s, str) or len(s) != length:
        return False
    return all(ch in NT_INDEX for ch in s)


def context_34_to_31(ctx: str) -> str:
    """Slice a 34-nt crisprware context to Chen 2025's 31-nt input window.

    Verified against `data/EnDeepCpf1/dataset/HEK_HCT_plasmid.csv` and
    `EnDeepCpf1_train.csv` in the capsule: the published `seq` column is
    exactly `seq34[1:32]`.
    """
    return ctx[1:32]


def one_hot_encode_31(seqs_31: Sequence[str]) -> np.ndarray:
    """One-hot encode 31-nt strings to channels-last shape (N, 31, 4) float32.

    Encoding order is A=0, T=1, C=2, G=3 (Chen's convention, NOT DeepCpf1's
    A=0, C=1, G=2, T=3 order).
    """
    n = len(seqs_31)
    arr = np.zeros((n, SEQ_LEN, 4), dtype=np.float32)
    for i, seq in enumerate(seqs_31):
        if not isinstance(seq, str) or len(seq) != SEQ_LEN:
            continue
        for j, ch in enumerate(seq):
            idx = NT_INDEX.get(ch)
            if idx is not None:
                arr[i, j, idx] = 1.0
    return arr


def predict(
    seqs: Sequence[str],
    model=None,
    batch_size: int = 1024,
    weights_path: Optional[str] = None,
    input_length: int = CONTEXT_LEN,
) -> List[float]:
    """Score sequences. Invalid (non-ACGT or wrong length) sequences return NaN.

    `input_length` defaults to 34 (crisprware standard Cas12a context). Pass
    `input_length=31` if you're feeding pre-sliced Chen-format strings directly
    (e.g. the test fixtures from the capsule's CSVs).
    """
    if model is None:
        model = load_model(weights_path)

    n = len(seqs)
    out = [float("nan")] * n
    valid_idx: List[int] = []
    seqs_31: List[str] = []

    for i, s in enumerate(seqs):
        if not is_valid_seq(s, length=input_length):
            continue
        s_31 = context_34_to_31(s) if input_length == CONTEXT_LEN else s
        if is_valid_seq(s_31, length=SEQ_LEN):
            valid_idx.append(i)
            seqs_31.append(s_31)

    if not valid_idx:
        return out

    x = one_hot_encode_31(seqs_31)
    y = model.predict(x, batch_size=min(batch_size, len(seqs_31)), verbose=0).ravel()
    for slot, v in zip(valid_idx, y):
        out[slot] = float(v)
    return out


def compute_enseq_deepcpf1_scores(
    sequences: Iterable[str],
    threads: int = 1,
    chunk_size: int = 100_000,
    weights_path: Optional[str] = None,
) -> List[float]:
    """Chunked entry mirroring `compute_rs3_scores` / `compute_deepcpf1_scores`."""
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
    out_col: str = "enseq_deepcpf1_score",
    batch_size: int = 1024,
) -> None:
    """Read a TSV, score `seq_col`, append `out_col`, write TSV."""
    import pandas as pd

    df = pd.read_csv(input_path, sep="\t", dtype=str)
    if df.columns[0].startswith("#"):
        df.rename(columns={df.columns[0]: df.columns[0].lstrip("#")}, inplace=True)
    if seq_col not in df.columns:
        raise ValueError(f"Column {seq_col!r} not found. Available: {list(df.columns)}")
    seqs = df[seq_col].fillna("").tolist()
    scores = predict(seqs, batch_size=batch_size)
    df[out_col] = [None if (s != s) else round(float(s), 8) for s in scores]
    df.to_csv(output_path, sep="\t", index=False)
