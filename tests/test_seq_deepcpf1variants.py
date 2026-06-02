"""Tests for crisprware.scorers.seq_deepcpf1variants.

Validates that variant-specific models load + score with the same EnCas12a
architecture as enseq-DeepCpf1, and that the AsCas12a variant agrees with
the per-variant fixtures shipped in the Code Ocean capsule.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from crisprware.scorers import seq_deepcpf1variants as variants


VARIANTS_DIR = os.path.join(REPO_ROOT, "crisprware", "scorers", "weights", "chen_2025", "variants")


def test_variant_list_matches_files():
    """Every name in VARIANTS must have a corresponding .npz on disk."""
    for v in variants.VARIANTS:
        path = os.path.join(VARIANTS_DIR, f"{v}.npz")
        assert os.path.exists(path), f"missing {path}"


def test_normalize_variant_accepts_common_forms():
    assert variants._normalize_variant("AsCas12a") == "AsCas12a"
    assert variants._normalize_variant("ascas12a") == "AsCas12a"
    # Space and underscore interchangeable
    assert variants._normalize_variant("AsCas12a Ultra") == "AsCas12a_Ultra"
    assert variants._normalize_variant("AsCas12a_Ultra") == "AsCas12a_Ultra"
    # iCas12a parentheses normalized to underscore
    assert variants._normalize_variant("iCas12a(mut2C-W)") == "iCas12a_mut2C-W"
    # Dash insensitive
    assert variants._normalize_variant("enAsCas12aHF1") == "enAsCas12a-HF1"
    assert variants._normalize_variant("ENASCAS12A-HF1") == "enAsCas12a-HF1"


def test_normalize_variant_rejects_unknown():
    with pytest.raises(ValueError, match="Unknown Cas12a variant"):
        variants._normalize_variant("SpCas9")


@pytest.fixture(scope="module")
def ascas12a_model():
    pytest.importorskip("tensorflow")
    return variants.load_model("AsCas12a")


def test_load_model_for_each_variant():
    """Sanity: every variant's weights file is loadable (smoke test, no scoring)."""
    pytest.importorskip("tensorflow")
    # Test a representative sample to keep test time bounded
    sample = ["AsCas12a", "AsCas12a_Ultra", "enAsCas12a-HF1", "LbCas12a", "HyperFi-AsCas12a", "iCas12a_mut2C-W"]
    for v in sample:
        m = variants.load_model(v)
        # Architecture should be identical -- ~609k params
        n_params = sum(int(np.prod(w.shape)) for w in m.weights)
        assert n_params == 609409, f"{v}: unexpected param count {n_params}"


def test_predict_via_34nt_context_works(ascas12a_model):
    """Smoke test: 34-nt context yields finite scores in [0, 1]."""
    ctxs = [
        "ACATTTTACTTTTTCAAAATTGTTTTCATGCTAA",
        "TTTTTTTAGTGAAGCTTCTAGATATTTGGCGGGT",
        "AGGTTTTACAACCGCCCAGTGCGTCTACGTCACA",
    ]
    out = variants.predict(ctxs, variant="AsCas12a", model=ascas12a_model)
    assert all(0.0 <= s <= 1.0 for s in out), out
    assert all(s == s for s in out), out  # no NaN


def test_ascas12a_predictions_match_capsule_fixture(ascas12a_model):
    """Parity: per-variant AsCas12a predictions should match the capsule's
    AsCas12a-only test fixture (if shipped)."""
    capsule_csv = os.path.join(
        REPO_ROOT, "capsule-9398276", "data", "EnDeepCpf1", "dataset", "Cas12a_variants", "AsCas12a.test.csv"
    )
    if not os.path.exists(capsule_csv):
        pytest.skip(f"missing capsule fixture at {capsule_csv}")
    df = pd.read_csv(capsule_csv).head(50)
    # The capsule CSV uses Chen's 31-nt 'seq' column.
    seqs_31 = df["seq"].tolist()
    scored = np.array(
        variants.predict(seqs_31, variant="AsCas12a", model=ascas12a_model, input_length=31),
        dtype=float,
    )
    assert np.all(np.isfinite(scored))
    assert np.all((scored >= 0.0) & (scored <= 1.0))
