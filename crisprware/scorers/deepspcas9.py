"""DeepSpCas9 on-target activity scorer (Kim et al. 2019, Sci Adv).

TF2/Keras port of the original Python 2.7 + TF 1.4 inception-style CNN
from Kim, Kim, Kim, Lee 2019 (DOI 10.1126/sciadv.aax9249, source MIT).
Weights extracted once from the upstream TF1 checkpoint via
`tools/extract_deepspcas9_weights.py` and bundled as a numpy .npz so the
runtime depends only on TF2.

Architecture
------------
Input:   one-hot (1, 30, 4) NHWC, A=0, C=1, G=2, T=3.
Conv:    3 parallel Conv2D inception branches with kernel sizes 3/5/7
         and filter counts 100/70/40, followed by AveragePooling2D
         (1, 2) stride 2 with 'same' padding.
Concat:  the three pooled outputs flatten and concatenate to a 2790-dim
         vector (1400 + 910 + 480 from kernels 3, 5, 7).
FC:      Dense(80) -> Dense(60) -> Dense(1, linear).

The 30-nt window is 4 nt 5' flank + 20 nt protospacer + 3 nt PAM (NGG)
+ 3 nt 3' flank, identical to crisprware's existing SpCas9 context
emitted by `generate_guides` with `--context_window 4 3 --pam NGG
--sgRNA_length 20` (no `--pam_5_prime`).

Output is an unbounded regression value (~[0, 100]) representing
predicted on-target activity. NOT normalized to [0, 1].
"""
from __future__ import annotations

import os
from importlib import resources
from typing import Iterable, List, Optional, Sequence

import numpy as np

SEQ_LEN = 30
NT_INDEX = {
    "A": 0, "C": 1, "G": 2, "T": 3,
    "a": 0, "c": 1, "g": 2, "t": 3,
}

_DEFAULT_WEIGHTS = os.path.join("kim_2019_deepspcas9", "deepspcas9.npz")

_LAYER_KEYS = [
    ("conv1", "conv1_w", "conv1_b"),
    ("conv2", "conv2_w", "conv2_b"),
    ("conv3", "conv3_w", "conv3_b"),
    ("fc1",   "fc1_w",   "fc1_b"),
    ("fc2",   "fc2_w",   "fc2_b"),
    ("out",   "out_w",   "out_b"),
]


def _weights_path(name: str = _DEFAULT_WEIGHTS) -> str:
    try:
        return str(resources.files(__package__).joinpath("weights", name))
    except (AttributeError, ModuleNotFoundError):
        return os.path.join(os.path.dirname(__file__), "weights", name)


def _build_model():
    """Build the TF2 Keras inception-CNN matching Kim 2019's DeepCas9 class."""
    import tensorflow as tf
    L = tf.keras.layers

    inp = L.Input(shape=(1, SEQ_LEN, 4), name="input")
    branches = []
    for kernel, n_filters, conv_name in [(3, 100, "conv1"), (5, 70, "conv2"), (7, 40, "conv3")]:
        x = L.Conv2D(
            n_filters, (1, kernel),
            padding="valid", activation="relu", name=conv_name,
        )(inp)
        x = L.AveragePooling2D(
            (1, 2), strides=(1, 2), padding="same", name=f"pool_{conv_name[-1]}"
        )(x)
        x = L.Flatten(name=f"flat_{conv_name[-1]}")(x)
        branches.append(x)
    concat = L.Concatenate(axis=-1, name="concat")(branches)
    x = L.Dense(80, activation="relu", name="fc1")(concat)
    x = L.Dense(60, activation="relu", name="fc2")(x)
    out = L.Dense(1, activation="linear", name="out")(x)
    return tf.keras.Model(inp, out, name="deepspcas9")


def load_model(weights_path: Optional[str] = None):
    """Build the model and load weights from a numpy .npz."""
    model = _build_model()
    path = weights_path or _weights_path()
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"DeepSpCas9 weights not found at {path}. "
            f"Run `tools/extract_deepspcas9_weights.py` once to populate it."
        )
    w = np.load(path)
    for layer_name, w_key, b_key in _LAYER_KEYS:
        model.get_layer(layer_name).set_weights([w[w_key], w[b_key]])
    return model


def is_valid_seq(s: object, length: int = SEQ_LEN) -> bool:
    if not isinstance(s, str) or len(s) != length:
        return False
    return all(ch in NT_INDEX for ch in s)


def one_hot_encode(seqs: Sequence[str]) -> np.ndarray:
    """One-hot encode 30-nt sequences to (N, 1, 30, 4) float32.

    Encoding follows Kim 2019 verbatim: A=0, C=1, G=2, T=3.
    """
    n = len(seqs)
    arr = np.zeros((n, 1, SEQ_LEN, 4), dtype=np.float32)
    for i, seq in enumerate(seqs):
        if not isinstance(seq, str) or len(seq) != SEQ_LEN:
            continue
        for j, ch in enumerate(seq):
            idx = NT_INDEX.get(ch)
            if idx is not None:
                arr[i, 0, j, idx] = 1.0
    return arr


def predict(
    seqs: Sequence[str],
    model=None,
    batch_size: int = 1024,
    weights_path: Optional[str] = None,
) -> List[float]:
    """Score sequences. Invalid (non-30-nt-ACGT) sequences return NaN."""
    if model is None:
        model = load_model(weights_path)

    n = len(seqs)
    out = [float("nan")] * n
    valid_idx = [i for i, s in enumerate(seqs) if is_valid_seq(s)]
    if not valid_idx:
        return out

    valid_seqs = [seqs[i] for i in valid_idx]
    x = one_hot_encode(valid_seqs)
    y = model.predict(x, batch_size=min(batch_size, len(valid_seqs)), verbose=0).ravel()
    for slot, v in zip(valid_idx, y):
        out[slot] = float(v)
    return out


def compute_deepspcas9_scores(
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
        scores.extend(predict(seqs[i:i + chunk_size], model=model))
    return scores


def score_file(
    input_path: str,
    output_path: str,
    seq_col: str = "context",
    out_col: str = "deepspcas9_score",
    batch_size: int = 1024,
) -> None:
    """Read a TSV, score `seq_col`, append `out_col`, write TSV."""
    import pandas as pd

    df = pd.read_csv(input_path, sep="\t", dtype=str)
    if df.columns[0].startswith("#"):
        df.rename(columns={df.columns[0]: df.columns[0].lstrip("#")}, inplace=True)
    if seq_col not in df.columns:
        raise ValueError(f"Column {seq_col!r} not found.")
    seqs = df[seq_col].fillna("").tolist()
    scores = predict(seqs, batch_size=batch_size)
    df[out_col] = [None if (s != s) else round(float(s), 8) for s in scores]
    df.to_csv(output_path, sep="\t", index=False)
