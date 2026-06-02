"""Parity tests for crisprware.scorers.deepspcas9.

Validates the TF2 port against ground-truth scores from the original
Python 2.7 + TF 1.4 reference in crisprScore's basilisk env. Reference
output (`deepspcas9_ref_output.txt`) was generated once via
`getDeepSpCas9Scores.py` on the 7-sequence fixture `deepspcas9_test_input.tsv`.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from crisprware.scorers import deepspcas9

WEIGHTS = os.path.join(REPO_ROOT, "crisprware", "scorers", "weights", "kim_2019_deepspcas9", "deepspcas9.npz")
TEST_INPUT = os.path.join(REPO_ROOT, "tests", "test_data", "deepspcas9_test_input.tsv")
REF_OUTPUT = os.path.join(REPO_ROOT, "tests", "test_data", "deepspcas9_ref_output.txt")

weights_required = pytest.mark.skipif(
    not os.path.exists(WEIGHTS),
    reason=f"missing {WEIGHTS} -- run `tools/extract_deepspcas9_weights.py` first",
)


def test_is_valid_seq():
    assert deepspcas9.is_valid_seq("A" * 30)
    assert not deepspcas9.is_valid_seq("A" * 29)
    assert not deepspcas9.is_valid_seq("N" + "A" * 29)


def test_one_hot_encode_shape_and_layout():
    arr = deepspcas9.one_hot_encode(["ACGT" * 7 + "AC"])
    assert arr.shape == (1, 1, 30, 4)
    assert arr.dtype == np.float32
    # A=0, C=1, G=2, T=3 (NOT Chen's ATCG order)
    np.testing.assert_array_equal(arr[0, 0, 0, :], [1, 0, 0, 0])
    np.testing.assert_array_equal(arr[0, 0, 1, :], [0, 1, 0, 0])
    np.testing.assert_array_equal(arr[0, 0, 2, :], [0, 0, 1, 0])
    np.testing.assert_array_equal(arr[0, 0, 3, :], [0, 0, 0, 1])


@pytest.fixture(scope="module")
def model():
    pytest.importorskip("tensorflow")
    if not os.path.exists(WEIGHTS):
        pytest.skip(f"missing {WEIGHTS}")
    return deepspcas9.load_model(WEIGHTS)


@weights_required
def test_model_param_count(model):
    """Architecture should produce exactly 240k params."""
    n_params = sum(int(np.prod(w.shape)) for w in model.weights)
    # 100*(1*3*4) + 100 + 70*(1*5*4) + 70 + 40*(1*7*4) + 40 +
    # 2790*80 + 80 + 80*60 + 60 + 60*1 + 1 = 232,131
    assert n_params == 232131, f"unexpected param count {n_params}"


@weights_required
def test_predict_matches_tf1_reference(model):
    """All 7 fixture sequences must match the TF1 reference within float32 noise."""
    import pandas as pd

    df = pd.read_csv(TEST_INPUT, sep="\t")
    seqs = df["sequence"].tolist()
    ref = np.loadtxt(REF_OUTPUT, dtype=float)
    assert len(seqs) == len(ref), (len(seqs), len(ref))

    scored = np.array(deepspcas9.predict(seqs, model=model), dtype=float)
    max_diff = float(np.max(np.abs(scored - ref)))
    assert max_diff < 1e-3, (
        f"max abs diff {max_diff} exceeds float32 noise tolerance.\nscored={scored.tolist()}\nref={ref.tolist()}"
    )


@weights_required
def test_predict_invalid_seq_returns_nan(model):
    out = deepspcas9.predict(["A" * 29, "A" * 30], model=model)
    assert np.isnan(out[0])
    assert not np.isnan(out[1])
