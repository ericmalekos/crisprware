"""Unit tests for crisprware.scorers.deepcpf1."""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from crisprware.scorers import deepcpf1

GOLDEN_TSV = os.path.join(
    REPO_ROOT, "tests", "test_data", "cas12a_scores", "deepcpf1_enpam_sgRNA.1.tsv"
)
WEIGHTS = os.path.join(
    REPO_ROOT, "crisprware", "scorers", "weights", "Seq_deepCpf1_weights.h5"
)


def test_is_valid_seq():
    assert deepcpf1.is_valid_seq("A" * 34)
    assert deepcpf1.is_valid_seq("acgt" * 8 + "ac")
    assert not deepcpf1.is_valid_seq("A" * 33)
    assert not deepcpf1.is_valid_seq("A" * 35)
    assert not deepcpf1.is_valid_seq("N" + "A" * 33)
    assert not deepcpf1.is_valid_seq(None)
    assert not deepcpf1.is_valid_seq(34)


def test_one_hot_encode_shape_and_layout():
    arr = deepcpf1.one_hot_encode(["ACGT" * 8 + "AC"])
    assert arr.shape == (1, 34, 4)
    assert arr.dtype == np.int8
    # First four positions should be A, C, G, T
    np.testing.assert_array_equal(arr[0, 0], [1, 0, 0, 0])
    np.testing.assert_array_equal(arr[0, 1], [0, 1, 0, 0])
    np.testing.assert_array_equal(arr[0, 2], [0, 0, 1, 0])
    np.testing.assert_array_equal(arr[0, 3], [0, 0, 0, 1])


def test_one_hot_encode_lowercase_and_invalid_chars():
    arr = deepcpf1.one_hot_encode(["aN" + "A" * 32])
    assert arr.shape == (1, 34, 4)
    np.testing.assert_array_equal(arr[0, 0], [1, 0, 0, 0])
    # N -> all zeros
    np.testing.assert_array_equal(arr[0, 1], [0, 0, 0, 0])


@pytest.fixture(scope="module")
def model():
    pytest.importorskip("tensorflow")
    assert os.path.exists(WEIGHTS), f"Missing weights at {WEIGHTS}"
    return deepcpf1.load_model(WEIGHTS)


def test_load_legacy_weights_shapes(model):
    """Conv1D kernel must come out as (5, 4, 80); dense layers as documented."""
    # Find the Conv1D layer
    conv_layers = [l for l in model.layers if l.__class__.__name__ == "Conv1D"]
    assert len(conv_layers) == 1
    assert tuple(conv_layers[0].weights[0].shape) == (5, 4, 80)
    # Dense layers in order: 1200->80, 80->40, 40->40, 40->1
    dense_shapes = [
        tuple(l.weights[0].shape)
        for l in model.layers
        if l.__class__.__name__ == "Dense"
    ]
    assert dense_shapes == [(1200, 80), (80, 40), (40, 40), (40, 1)]


def test_predict_invalid_seq_returns_nan(model):
    out = deepcpf1.predict(["A" * 33, "A" * 34], model=model)
    assert np.isnan(out[0])
    assert not np.isnan(out[1])


def test_predict_matches_golden_tsv(model):
    """Top 10 rows of the user's golden fixture must score within tolerance."""
    df = pd.read_csv(GOLDEN_TSV, sep="\t").head(10)
    seqs = df["context"].tolist()
    golden = df["deepcpf1_score"].astype(float).tolist()
    scored = deepcpf1.predict(seqs, model=model)
    diffs = [abs(s - g) for s, g in zip(scored, golden)]
    max_diff = max(diffs)
    assert max_diff < 1e-2, (
        f"DeepCpf1 port diverges from golden by {max_diff} on top-10 fixture. "
        f"Scored={scored}, golden={golden}"
    )


def test_predict_spearman_high_correlation(model):
    """Across the full first-chunk fixture, Spearman vs golden must be ~1."""
    from scipy.stats import spearmanr

    df = pd.read_csv(GOLDEN_TSV, sep="\t")
    seqs = df["context"].tolist()
    golden = df["deepcpf1_score"].astype(float).to_numpy()
    scored = np.array(deepcpf1.predict(seqs, model=model))
    rho, _ = spearmanr(scored, golden)
    assert rho >= 0.9999, f"Spearman {rho} below 0.9999 vs golden"
    mean_abs_diff = float(np.mean(np.abs(scored - golden)))
    assert mean_abs_diff <= 1e-4, f"Mean abs diff {mean_abs_diff} > 1e-4"
