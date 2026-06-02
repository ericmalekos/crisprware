"""
DeepCpf1 (Seq-deepCpf1) on-target scorer, TF2 port.

Architecture from Kim et al., Nat Biotechnol 2018 (PMID 29431740):
    Input (34 x 4) -> Conv1D(80, k=5, relu) -> AvgPool1D(2) ->
    Flatten -> Dropout(0.3) -> Dense(80, relu) -> Dropout(0.3) ->
    Dense(40, relu) -> Dropout(0.3) -> Dense(40, relu) -> Dropout(0.3) ->
    Dense(1, linear)

Loads weights either from the user's legacy Theano-Keras .h5
(`Seq_deepCpf1_weights.h5`, with the 4-D Conv1D kernel layout
`(kernel, 1, in_ch, out_ch)`) or a TF2-shaped .h5 (kernel `(k, in, out)`).
Layer order is taken from the h5 `layer_names` attribute when present;
otherwise from the h5 group keys.

Public surface mirrors `rs3` (see `score_guides.compute_rs3_scores`):
    one_hot_encode, predict, compute_deepcpf1_scores, score_file.
"""

from __future__ import annotations

import os
from importlib import resources
from typing import Iterable, List, Optional, Sequence

import numpy as np

SEQ_LEN = 34
NT_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}

_DEFAULT_WEIGHTS = "Seq_deepCpf1_weights.h5"


def _weights_path() -> str:
    """Locate the bundled weights file inside the installed package."""
    try:
        return str(resources.files(__package__).joinpath("weights", _DEFAULT_WEIGHTS))
    except (AttributeError, ModuleNotFoundError):
        return os.path.join(os.path.dirname(__file__), "weights", _DEFAULT_WEIGHTS)


def _build_model():
    """Build the Seq-deepCpf1 architecture in TF2 / Keras 2.x."""
    import tensorflow as tf

    L = tf.keras.layers

    inp = L.Input(shape=(SEQ_LEN, 4), name="input_1")
    x = L.Conv1D(80, 5, activation="relu", padding="valid", name="conv1d_1")(inp)
    x = L.AveragePooling1D(2, name="avg_pool_1")(x)
    x = L.Flatten(name="flatten_1")(x)
    x = L.Dropout(0.3, name="dropout_1")(x)
    x = L.Dense(80, activation="relu", name="dense_1")(x)
    x = L.Dropout(0.3, name="dropout_2")(x)
    x = L.Dense(40, activation="relu", name="dense_2")(x)
    x = L.Dropout(0.3, name="dropout_3")(x)
    x = L.Dense(40, activation="relu", name="dense_3")(x)
    x = L.Dropout(0.3, name="dropout_4")(x)
    out = L.Dense(1, activation="linear", name="dense_out")(x)
    return tf.keras.Model(inp, out, name="seq_deepcpf1")


def _legacy_layer_groups(f) -> List[str]:
    """Return h5 group names in layer order.

    Old Keras stored layer order under the file-level `layer_names` attribute;
    newer files keep it under `model_weights.layer_names`. If neither, fall
    back to insertion order of the top-level groups.
    """
    if "layer_names" in f.attrs:
        return [n.decode() if isinstance(n, bytes) else n for n in f.attrs["layer_names"]]
    if "model_weights" in f and "layer_names" in f["model_weights"].attrs:
        return [n.decode() if isinstance(n, bytes) else n for n in f["model_weights"].attrs["layer_names"]]
    return list(f.keys())


def _extract_weights_pair(f, group_name: str):
    """Return (kernel, bias) numpy arrays from an h5 layer group, or (None, None)."""
    g = f[group_name]

    candidates = list(g.keys())
    if "model_weights" in f and len(candidates) == 0:
        g = f["model_weights"][group_name]
        candidates = list(g.keys())

    inner = None
    for k in candidates:
        if hasattr(g[k], "shape"):
            inner = g
            break
        if k == group_name:
            inner = g[k]
            break
    if inner is None and candidates:
        inner = g[candidates[0]]

    kernel = bias = None
    for k in inner.keys() if hasattr(inner, "keys") else []:
        if hasattr(inner[k], "shape"):
            kl = k.lower()
            if kl.endswith("kernel:0") or kl.endswith("_w") or kl == "kernel":
                kernel = np.array(inner[k])
            elif kl.endswith("bias:0") or kl.endswith("_b") or kl == "bias":
                bias = np.array(inner[k])
    return kernel, bias


def _load_legacy_h5(model, path: str) -> None:
    """Iterate h5 groups in layer order and set weights into the TF2 model.

    Handles two h5 layouts for the Conv1D kernel:
      - `(k, 1, in_ch, out_ch)` from Theano-Keras 1.x (4-D Conv1D, stored in
        Theano's mathematical-convolution convention). Squeezed to
        `(k, in_ch, out_ch)` and flipped along axis 0 to match TF's
        cross-correlation Conv1D.
      - `(k, in_ch, out_ch)` from TF Keras — loaded verbatim.
    """
    import h5py

    target_layers = [layer for layer in model.layers if layer.weights]
    target_iter = iter(target_layers)

    with h5py.File(path, "r") as f:
        for group_name in _legacy_layer_groups(f):
            kernel, bias = _extract_weights_pair(f, group_name)
            if kernel is None and bias is None:
                continue

            try:
                tgt = next(target_iter)
            except StopIteration:
                break

            expected_shape = tgt.weights[0].shape
            theano_conv = False

            if kernel.ndim == 4 and kernel.shape[1] == 1:
                kernel = kernel.squeeze(axis=1)
                theano_conv = True
            elif kernel.ndim == 4 and kernel.shape[2] == 1:
                kernel = kernel.squeeze(axis=2)
                theano_conv = True

            if theano_conv and tgt.__class__.__name__ == "Conv1D":
                kernel = kernel[::-1].copy()

            if tuple(kernel.shape) != tuple(expected_shape):
                raise ValueError(
                    f"Weight shape mismatch for layer {tgt.name}: "
                    f"h5 group {group_name} kernel shape {kernel.shape} "
                    f"vs target {tuple(expected_shape)}"
                )

            tgt.set_weights([kernel, bias])


def _is_theano_keras1_h5(path: str) -> bool:
    """Return True if the h5 file has the Theano-Keras 1.x layout.

    Detected by either (a) a file-level `layer_names` attr (instead of the
    newer `model_weights/layer_names`), or (b) a Conv1D weight tensor with
    4-D shape `(kernel, 1, in_ch, out_ch)`.
    """
    import h5py

    with h5py.File(path, "r") as f:
        if "layer_names" in f.attrs:
            return True
        for k in f.keys():
            if "conv" in k.lower():
                grp = f[k]
                for name in grp.keys() if hasattr(grp, "keys") else []:
                    v = grp[name]
                    if hasattr(v, "shape") and v.ndim == 4:
                        return True
    return False


def load_model(weights_path: Optional[str] = None):
    """Build the model and load weights from `weights_path` (or the bundled file)."""
    model = _build_model()
    path = weights_path or _weights_path()
    if _is_theano_keras1_h5(path):
        _load_legacy_h5(model, path)
    else:
        model.load_weights(path)
    return model


def is_valid_seq(s: object) -> bool:
    """Return True if `s` is a 34-character ACGT string (case-insensitive)."""
    if not isinstance(s, str) or len(s) != SEQ_LEN:
        return False
    return all(ch in NT_INDEX for ch in s)


def one_hot_encode(seqs: Sequence[str]) -> np.ndarray:
    """One-hot encode 34-nt sequences to shape (N, 34, 4), int8.

    Encoding: A->[1,0,0,0], C->[0,1,0,0], G->[0,0,1,0], T->[0,0,0,1].
    Non-ACGT positions remain all-zero (treated as N).
    """
    n = len(seqs)
    arr = np.zeros((n, SEQ_LEN, 4), dtype=np.int8)
    for i, seq in enumerate(seqs):
        if not isinstance(seq, str):
            continue
        for j, ch in enumerate(seq[:SEQ_LEN]):
            idx = NT_INDEX.get(ch)
            if idx is not None:
                arr[i, j, idx] = 1
    return arr


def predict(
    seqs: Sequence[str],
    model=None,
    batch_size: int = 1024,
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

    valid_seqs = [seqs[i] for i in valid_idx]
    x = one_hot_encode(valid_seqs).astype(np.float32)
    y = model.predict(x, batch_size=min(batch_size, len(valid_seqs)), verbose=0).ravel()
    for slot, v in zip(valid_idx, y):
        out[slot] = float(v)
    return out


def compute_deepcpf1_scores(
    sequences: Iterable[str],
    threads: int = 1,
    chunk_size: int = 100_000,
    weights_path: Optional[str] = None,
) -> List[float]:
    """Chunked entry point mirroring `compute_rs3_scores` in score_guides.py.

    `threads` is accepted for signature symmetry with the rs3 path; TF
    handles parallelism internally via its thread pool. Set the env vars
    `TF_NUM_INTRAOP_THREADS` / `TF_NUM_INTEROP_THREADS` before importing TF
    if you need to constrain it.
    """
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
    out_col: str = "deepcpf1_score",
    batch_size: int = 1024,
) -> None:
    """Read a TSV, score `seq_col` with DeepCpf1, append `out_col`, write TSV.

    Drop-in for the legacy `parasol_scripts/DeepCpf1_seq.py` (same flags).
    """
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
