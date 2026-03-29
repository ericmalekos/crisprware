#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DeepCpf1 TSV scorer.

Reads a TSV file, scores sequences in a specified column (default: "context")
using the DeepCpf1 model, and writes an output TSV with a new score column.

Usage:
    python DeepCpf1_score.py \
      -i sgRNAs/sgRNA.6.tsv \
      -o sgRNAs/sgRNA.6.deepcpf1.tsv \
      --seq-col context \
      --out-col deepcpf1_score
"""

from __future__ import print_function
import os
import sys
import argparse

import numpy as np
import pandas as pd

# Keras/Theano (pinned in your env)
os.environ.setdefault("KERAS_BACKEND", "theano")  # must be before any `from keras...` imports
from keras import backend as K
from keras.models import Model
from keras.layers import Input
from keras.layers.core import Dense, Dropout, Flatten
from keras.layers.convolutional import Convolution1D, AveragePooling1D

# --- Config ---
SEQ_LEN = 34  # DeepCpf1 expects 34-nt context
dir_path = os.path.dirname(os.path.realpath(__file__))
weights_path = os.path.join(dir_path, "weights/Seq_deepCpf1_weights.h5")


# Compatibility: reload is builtin in py2; importlib.reload in py3
try:
    from importlib import reload  # py3
except Exception:
    pass  # py2 has builtin reload


def set_keras_backend(backend):
    """Force a specific Keras backend (e.g., 'theano')."""
    if K.backend() != backend:
        os.environ["KERAS_BACKEND"] = backend
        # Quietly reload the backend module
        _stderr = sys.stderr
        try:
            sys.stderr = open(os.devnull, "w")
            reload(K)
        finally:
            sys.stderr = _stderr
        assert K.backend() == backend, "Failed to switch Keras backend to {}".format(backend)


def build_deepcpf1_model():
    """Construct the DeepCpf1 model graph and load pretrained weights."""
    if K.backend() != "theano":
        print("ERROR: Keras backend is '{}', but DeepCpf1 requires 'theano'.".format(K.backend()), file=sys.stderr)
        print("       Set the backend or use the provided set_keras_backend('theano').", file=sys.stderr)
    inp = Input(shape=(SEQ_LEN, 4))
    x = Convolution1D(80, 5, activation="relu")(inp)
    x = AveragePooling1D(2)(x)
    x = Flatten()(x)
    x = Dropout(0.3)(x)
    x = Dense(80, activation="relu")(x)
    x = Dropout(0.3)(x)
    x = Dense(40, activation="relu")(x)
    x = Dropout(0.3)(x)
    x = Dense(40, activation="relu")(x)
    x = Dropout(0.3)(x)
    out = Dense(1, activation="linear")(x)
    model = Model(inputs=[inp], outputs=[out])
    model.load_weights(weights_path)
    return model


def preprocess_sequences(seq_list):
    """
    One-hot encode sequences to shape (N, 34, 4).

    Parameters
    ----------
    seq_list : list[str]
        List of 34-nt strings composed of A/C/G/T (case-insensitive).

    Returns
    -------
    numpy.ndarray
        One-hot matrix of shape (N, 34, 4), dtype int.
    """
    data_n = len(seq_list)
    SEQ = np.zeros((data_n, SEQ_LEN, 4), dtype=int)
    for l, seq in enumerate(seq_list):
        for i, ch in enumerate(seq):
            if ch in "Aa":
                SEQ[l, i, 0] = 1
            elif ch in "Cc":
                SEQ[l, i, 1] = 1
            elif ch in "Gg":
                SEQ[l, i, 2] = 1
            elif ch in "Tt":
                SEQ[l, i, 3] = 1
    return SEQ


def is_valid_seq(s):
    """Return True if s is a 34-nt string with only A/C/G/T (case-insensitive)."""
    if not isinstance(s, str) or len(s) != SEQ_LEN:
        return False
    for ch in s:
        if ch.upper() not in ("A", "C", "G", "T"):
            return False
    return True


def score_file(input_path, output_path, seq_col="context", out_col="deepcpf1_score", batch_size=1024):
    """
    Read TSV, score DeepCpf1 on seq_col, and write TSV with out_col.

    Parameters
    ----------
    input_path : str
        Path to input TSV (must have a header row).
    output_path : str
        Path to output TSV (overwritten if exists).
    seq_col : str
        Column containing sequences to score (default: "context").
    out_col : str
        Name of the new column for scores (default: "deepcpf1_score").
    batch_size : int
        Batch size for prediction (default: 1024).
    """
    df = pd.read_csv(input_path, sep="\t", dtype=str)

    # Normalize first header if it begins with '#'
    if len(df.columns) > 0 and isinstance(df.columns[0], str) and df.columns[0].startswith("#"):
        df.rename(columns={df.columns[0]: df.columns[0].lstrip("#")}, inplace=True)

    if seq_col not in df.columns:
        raise ValueError("Column '{}' not found. Available: {}".format(seq_col, list(df.columns)))

    seq_series = df[seq_col].fillna("")

    # Validate sequences
    valid_mask = seq_series.apply(is_valid_seq)
    invalid_n = (~valid_mask).sum()
    if invalid_n:
        print("[DeepCpf1] Warning: {} rows have invalid sequences (not 34-nt A/C/G/T) — setting score to NaN."
              .format(int(invalid_n)), file=sys.stderr)

    # Prepare model once
    set_keras_backend("theano")
    model = build_deepcpf1_model()

    # Predict in batches for valid sequences
    scores = np.full(len(df), np.nan, dtype=float)
    valid_idx = np.flatnonzero(valid_mask.values)              # integer positions of valid rows
    valid_seqs = seq_series.iloc[valid_idx].tolist()

    for i in range(0, len(valid_seqs), batch_size):
        j0, j1 = i, min(i + batch_size, len(valid_seqs))
        batch = valid_seqs[j0:j1]
        X = preprocess_sequences(batch)
        y = model.predict([X], batch_size=min(batch_size, len(batch)), verbose=0).ravel()
        scores[valid_idx[j0:j1]] = y  

    # Round to 8 decimal places
    df[out_col] = [None if np.isnan(x) else round(float(x), 8) for x in scores]

    df.to_csv(output_path, sep="\t", index=False)


def parse_args():
    p = argparse.ArgumentParser(description="Score sequences in a TSV with DeepCpf1 and append results.")
    p.add_argument("-i", "--input", required=True, help="Input TSV path")
    p.add_argument("-o", "--output", required=True, help="Output TSV path")
    p.add_argument("--seq-col", default="context", help="Column containing sequences to score (default: context)")
    p.add_argument("--out-col", default="deepcpf1_score", help="Name of the output score column")
    p.add_argument("--batch-size", type=int, default=1024, help="Batch size for prediction (default: 1024)")
    return p.parse_args()


def main():
    args = parse_args()
    score_file(args.input, args.output, args.seq_col, args.out_col, args.batch_size)
    print("Wrote DeepCpf1-scored TSV to {}".format(args.output))


if __name__ == "__main__":
    main()
