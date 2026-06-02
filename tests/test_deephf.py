"""Parity tests for crisprware.scorers.deephf.

Validates the TF2/Keras port (loading the .hdf5 files preserved by
Bioconductor's crisprScoreData / ExperimentHub) against ground-truth
scores captured from the original Keras 2.1.6 + Python 3.6 reference
running in the basilisk env. The 11-dim bio-feature vector is
reproduced with vendored Python (ViennaRNA for stem/dG, Biopython for
melting temps, pure-Python lookup table for dG_binding).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from crisprware.scorers import deephf

WEIGHTS_DIR = os.path.join(
    REPO_ROOT, "crisprware", "scorers", "weights", "wang_2019_deephf"
)
REF_TSV = os.path.join(REPO_ROOT, "tests", "test_data", "deephf_ref_scores.tsv")


def _weights_present(variant: str) -> bool:
    fn = deephf._VARIANT_TO_WEIGHTS.get(variant)
    return fn is not None and os.path.exists(os.path.join(WEIGHTS_DIR, fn))


weights_required = pytest.mark.skipif(
    not _weights_present("wt_u6"),
    reason="DeepHF wt_u6 weights missing",
)


def test_is_valid_seq():
    assert deephf.is_valid_seq("A" * 23)
    assert not deephf.is_valid_seq("A" * 22)
    assert not deephf.is_valid_seq("N" + "A" * 22)


def test_encode_sequence_layout():
    """make_data: START token at position 0, A=2, T=3, C=4, G=5."""
    arr = deephf.encode_sequence(["ATCGATGCTGATGCTAGATAAGG"])
    assert arr.shape == (1, 22)
    assert arr.dtype == np.int32
    # Position 0 is START
    assert arr[0, 0] == 1
    # Position 1 onwards: ATCGATGCTGATGCTAGATAA (first 21 chars of input)
    expected = [2, 3, 4, 5, 2, 3, 5, 4, 3, 5, 2, 3, 5, 4, 3, 2, 5, 2, 3, 2, 2]
    np.testing.assert_array_equal(arr[0, 1:22], expected)


def test_dG_binding_lookup_table():
    """Pure dinucleotide free-energy lookup, matches feature_util.py:416."""
    # dG(aa) = -0.2 (only pair) + init 3.1 = 2.9
    assert deephf._dG_binding("aa") == pytest.approx(2.9, abs=1e-6)
    # All gc steps: ggcc gives 3.1 + dG(gg) + dG(gc) + dG(cc)
    # = 3.1 + -2.1 + -2.7 + -2.9 = -4.6
    assert deephf._dG_binding("ggcc") == pytest.approx(-4.6, abs=1e-6)


def test_gc_features():
    above, below, cnt = deephf._gc_features("A" * 21)
    assert cnt == 0
    assert above == 0
    assert below == 1


@pytest.fixture(scope="module", params=["wt_u6", "esp", "hf"])
def variant_model(request):
    pytest.importorskip("tensorflow")
    pytest.importorskip("RNA")
    if not _weights_present(request.param):
        pytest.skip(f"missing weights for {request.param}")
    return request.param, deephf.load_model(request.param)


def test_bio_features_match_reference():
    """compute_bio_features should produce the same 11-vector basilisk produces."""
    pytest.importorskip("RNA")
    if not _weights_present("wt_u6"):
        pytest.skip("weights not present")
    # Captured from basilisk on the same sequences:
    #   ATCGATGCTGATGCTAGATAAGG -> [0, -2.4, -21.8, -11.9, 0, 1, 8, 48.60039, -33.41291, 14.34474, -54.56225]
    #   GGAAGTCTGGAGTCTCCAGGTGG -> [0, -5.5, -26.2, -16.3, 1, 0, 12, 54.56942, -11.30664, 16.07473, -61.37201]
    feat = deephf.compute_bio_features([
        "ATCGATGCTGATGCTAGATAAGG",
        "GGAAGTCTGGAGTCTCCAGGTGG",
    ])
    expected = np.array([
        [0., -2.4, -21.8, -11.9, 0., 1., 8., 48.60038976, -33.41290629, 14.34473886, -54.56225157],
        [0., -5.5, -26.2, -16.3, 1., 0., 12., 54.56942238, -11.30663628, 16.07473406, -61.37201071],
    ], dtype=np.float32)
    np.testing.assert_allclose(feat, expected, atol=1e-3)


@weights_required
def test_predict_matches_basilisk_reference(variant_model):
    """All 3 variants must match the captured basilisk-env ground-truth."""
    variant, model = variant_model
    df = pd.read_csv(REF_TSV, sep="\t")
    seqs = df["sequence"].tolist()
    ref = df[variant].astype(float).to_numpy()
    scored = np.array(deephf.predict(seqs, variant=variant, model=model), dtype=float)
    max_diff = float(np.max(np.abs(scored - ref)))
    assert max_diff < 1e-4, (
        f"DeepHF {variant} diverges from basilisk reference by {max_diff}.\n"
        f"scored={scored.tolist()}, ref={ref.tolist()}"
    )


@weights_required
def test_predict_invalid_seq_returns_nan(variant_model):
    _, model = variant_model
    out = deephf.predict(["A" * 22, "A" * 23], model=model)
    assert np.isnan(out[0])
    assert not np.isnan(out[1])
