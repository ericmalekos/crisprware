"""Integration tests for the Cas12a on-target scoring branch of score_guides.main().

Drives the full pipeline against the golden Cas12a fixture (S. cerevisiae chrI
guides with pre-computed enpam_gb / deepcpf1 scores) and asserts the new flag
wiring produces matching output.
"""
from __future__ import annotations

import os
import sys
import tempfile
from argparse import Namespace

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

GOLDEN_TSV = os.path.join(
    REPO_ROOT, "tests", "test_data", "cas12a_scores", "deepcpf1_enpam_sgRNA.1.tsv"
)
ENPAM_GB_WEIGHTS = os.path.join(
    REPO_ROOT, "crisprware", "scorers", "weights", "enpam_gb.joblib"
)


def _build_input_bed(tmpdir: str, n_rows: int = 5) -> str:
    """Slice the golden TSV's first n rows into a fresh BED (without score cols)."""
    df = pd.read_csv(GOLDEN_TSV, sep="\t").head(n_rows)
    in_path = os.path.join(tmpdir, "cas12a_test.bed")
    cols = ["chr", "start", "stop", "id,sequence,pam,chromosome,position,sense", "context", "strand"]
    df[cols].to_csv(in_path, sep="\t", index=False)
    return in_path


def _make_args(grna_bed: str, output_directory: str, cas12a_scorer: str,
               cas12a_variant=None) -> Namespace:
    """Build the Namespace score_guides.main() expects."""
    return Namespace(
        grna_bed=grna_bed,
        output_directory=output_directory,
        threads=1,
        chunk_size=100,
        drop_duplicates=True,
        skip_rs3=True,
        skip_gs2=True,
        guidescan2_indices=[],
        tracr=None,
        alt_pams=None,
        mismatches=3,
        rna_bulges=0,
        dna_bulges=0,
        threshold=2,
        mode="succinct",
        keep_tmp=False,
        min_rs3=float("-inf"),
        cas12a_scorer=cas12a_scorer,
        cas12a_variant=cas12a_variant,
        min_deepcpf1=float("-inf"),
        min_enpam_gb=float("-inf"),
        min_enseq_deepcpf1=float("-inf"),
    )


def _read_scored_output(output_directory: str) -> pd.DataFrame:
    """Locate the score_guides output TSV (a .bed file actually) under output_directory."""
    matches = []
    for root, _, files in os.walk(output_directory):
        for fn in files:
            if fn.endswith(".bed") and "scoredgRNA" in fn:
                matches.append(os.path.join(root, fn))
    assert len(matches) == 1, f"expected one scoredgRNA bed, got {matches}"
    return pd.read_csv(matches[0], sep="\t")


def test_deepcpf1_branch_wired():
    """`--cas12a_scorer deepcpf1` adds a deepcpf1_score column and matches golden."""
    pytest.importorskip("tensorflow")
    from crisprware import score_guides

    with tempfile.TemporaryDirectory() as tmp:
        in_bed = _build_input_bed(tmp, n_rows=5)
        out_dir = os.path.join(tmp, "out")
        os.makedirs(out_dir, exist_ok=True)

        args = _make_args(in_bed, out_dir, cas12a_scorer="deepcpf1")
        score_guides.main(args)
        df = _read_scored_output(out_dir)

    assert "deepcpf1_score" in df.columns, df.columns.tolist()
    golden = pd.read_csv(GOLDEN_TSV, sep="\t").head(5)
    np.testing.assert_allclose(
        df["deepcpf1_score"].to_numpy(),
        golden["deepcpf1_score"].to_numpy(),
        atol=1e-2,
        err_msg="deepcpf1_score mismatch vs golden",
    )


@pytest.mark.skipif(
    not os.path.exists(ENPAM_GB_WEIGHTS),
    reason="missing enpam_gb.joblib — run `python tools/extract_enpam_gb.py` first",
)
def test_enpam_gb_branch_wired():
    """`--cas12a_scorer enpam_gb` adds an enpam_gb_score column and matches golden."""
    from crisprware import score_guides

    with tempfile.TemporaryDirectory() as tmp:
        in_bed = _build_input_bed(tmp, n_rows=5)
        out_dir = os.path.join(tmp, "out")
        os.makedirs(out_dir, exist_ok=True)

        args = _make_args(in_bed, out_dir, cas12a_scorer="enpam_gb")
        score_guides.main(args)
        df = _read_scored_output(out_dir)

    assert "enpam_gb_score" in df.columns, df.columns.tolist()
    golden = pd.read_csv(GOLDEN_TSV, sep="\t").head(5)
    np.testing.assert_allclose(
        df["enpam_gb_score"].to_numpy(),
        golden["enpam_gb_score"].to_numpy(),
        atol=1e-4,
    )


@pytest.mark.skipif(
    not os.path.exists(ENPAM_GB_WEIGHTS),
    reason="missing enpam_gb.joblib — run `python tools/extract_enpam_gb.py` first",
)
def test_both_branches_wired():
    """`--cas12a_scorer both` adds both score columns."""
    pytest.importorskip("tensorflow")
    from crisprware import score_guides

    with tempfile.TemporaryDirectory() as tmp:
        in_bed = _build_input_bed(tmp, n_rows=5)
        out_dir = os.path.join(tmp, "out")
        os.makedirs(out_dir, exist_ok=True)

        args = _make_args(in_bed, out_dir, cas12a_scorer="both")
        score_guides.main(args)
        df = _read_scored_output(out_dir)

    assert "deepcpf1_score" in df.columns
    assert "enpam_gb_score" in df.columns


def test_mutex_tracr_and_cas12a_scorer():
    """Setting --tracr and --cas12a_scorer together must raise."""
    from crisprware import score_guides

    args = _make_args(grna_bed="dummy", output_directory="/tmp", cas12a_scorer="enpam_gb")
    args.tracr = "Chen2013"
    args.skip_rs3 = False
    with pytest.raises(ValueError, match="mutually exclusive"):
        score_guides.main(args)


def test_enseq_deepcpf1_branch_wired():
    """`--cas12a_scorer enseq_deepcpf1` adds an enseq_deepcpf1_score column."""
    pytest.importorskip("tensorflow")
    from crisprware import score_guides

    with tempfile.TemporaryDirectory() as tmp:
        in_bed = _build_input_bed(tmp, n_rows=5)
        out_dir = os.path.join(tmp, "out")
        os.makedirs(out_dir, exist_ok=True)

        args = _make_args(in_bed, out_dir, cas12a_scorer="enseq_deepcpf1")
        score_guides.main(args)
        df = _read_scored_output(out_dir)

    assert "enseq_deepcpf1_score" in df.columns, df.columns.tolist()
    scores = df["enseq_deepcpf1_score"].astype(float).to_numpy()
    assert np.all(np.isfinite(scores))
    assert np.all((scores >= 0.0) & (scores <= 1.0))


def test_seq_deepcpf1variants_requires_variant():
    """`seq_deepcpf1variants` without --cas12a_variant must raise."""
    from crisprware import score_guides

    args = _make_args(grna_bed="dummy", output_directory="/tmp",
                      cas12a_scorer="seq_deepcpf1variants")
    with pytest.raises(ValueError, match="--cas12a_variant is required"):
        score_guides.main(args)


def test_seq_deepcpf1variants_AsCas12a_Ultra_branch_wired():
    """`--cas12a_scorer seq_deepcpf1variants --cas12a_variant AsCas12a_Ultra`
    adds a variant-specific score column."""
    pytest.importorskip("tensorflow")
    from crisprware import score_guides

    with tempfile.TemporaryDirectory() as tmp:
        in_bed = _build_input_bed(tmp, n_rows=5)
        out_dir = os.path.join(tmp, "out")
        os.makedirs(out_dir, exist_ok=True)

        args = _make_args(in_bed, out_dir,
                          cas12a_scorer="seq_deepcpf1variants",
                          cas12a_variant="AsCas12a_Ultra")
        score_guides.main(args)
        df = _read_scored_output(out_dir)

    expected_col = "seq_deepcpf1variants_AsCas12a_Ultra_score"
    assert expected_col in df.columns, df.columns.tolist()
