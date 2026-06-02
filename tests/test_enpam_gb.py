"""Unit tests for crisprware.scorers.enpam_gb.

Tests requiring the joblib model auto-skip if `crisprware/scorers/weights/enpam_gb.joblib`
hasn't been extracted yet. Run `python tools/extract_enpam_gb.py` once to populate it.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from crisprware.scorers import enpam_gb

GOLDEN_TSV = os.path.join(REPO_ROOT, "tests", "test_data", "cas12a_scores", "deepcpf1_enpam_sgRNA.1.tsv")
WEIGHTS = os.path.join(REPO_ROOT, "crisprware", "scorers", "weights", "enpam_gb.joblib")

weights_required = pytest.mark.skipif(
    not os.path.exists(WEIGHTS),
    reason=f"missing {WEIGHTS} — run `python tools/extract_enpam_gb.py` first",
)


def test_is_valid_seq():
    assert enpam_gb.is_valid_seq("A" * 34)
    assert enpam_gb.is_valid_seq("acgt" * 8 + "ac")
    assert not enpam_gb.is_valid_seq("A" * 33)
    assert not enpam_gb.is_valid_seq("N" + "A" * 33)


def test_featurize_shape_and_keys():
    df = enpam_gb.featurize(["A" * 34])
    expected_n = (
        1  # GC content
        + len(enpam_gb.NTS)  # Pos. Ind. 1mer  = 4
        + len(enpam_gb.NTS) ** 2  # Pos. Ind. 2mer  = 16
        + enpam_gb.CONTEXT_LEN * len(enpam_gb.NTS)  # Pos. Dep. 1mer = 34*4 = 136
        + (enpam_gb.CONTEXT_LEN - 1) * len(enpam_gb.NTS) ** 2  # Pos. Dep. 2mer = 33*16 = 528
        + 4  # Tm
    )
    assert df.shape == (1, expected_n), (df.shape, expected_n)
    assert "GC content" in df.columns
    assert "Tm, context" in df.columns


def test_guide_extraction():
    # PAM is positions 5-8 (1-indexed), guide is positions 9-31, downstream 32-34
    ctx = "AAAATTTV" + "C" * 23 + "GGG"  # 8 + 23 + 3 = 34
    assert len(ctx) == 34
    guide = enpam_gb._get_guide_sequence(ctx)
    assert guide == "C" * 23
    assert len(guide) == 23


@weights_required
def test_predict_matches_golden_top10():
    """Top 10 rows of the user's golden fixture must score within tolerance."""
    df = pd.read_csv(GOLDEN_TSV, sep="\t").head(10)
    seqs = df["context"].tolist()
    golden = df["enpam_gb_score"].astype(float).tolist()
    scored = enpam_gb.predict(seqs)
    diffs = [abs(s - g) for s, g in zip(scored, golden)]
    max_diff = max(diffs)
    assert max_diff < 1e-4, f"enPAM+GB port diverges from golden by {max_diff}. Scored={scored}, golden={golden}"


@weights_required
def test_predict_spearman_high_correlation():
    from scipy.stats import spearmanr

    df = pd.read_csv(GOLDEN_TSV, sep="\t")
    seqs = df["context"].tolist()
    golden = df["enpam_gb_score"].astype(float).to_numpy()
    scored = np.array(enpam_gb.predict(seqs))
    rho, _ = spearmanr(scored, golden)
    assert rho >= 0.9999, f"Spearman {rho} below 0.9999 vs golden"
    mean_abs_diff = float(np.mean(np.abs(scored - golden)))
    assert mean_abs_diff <= 1e-4, f"Mean abs diff {mean_abs_diff} > 1e-4"


@weights_required
def test_predict_nan_for_invalid():
    out = enpam_gb.predict(["A" * 33, "A" * 34])
    assert np.isnan(out[0])
    assert not np.isnan(out[1])
