#!/usr/bin/env python3
"""Baseline benchmark: score Kim 2018's endogenous test sheets with the
vendored DeepCpf1 and enPAM+GB ports, report Spearman / Pearson against
measured indel frequencies.

These numbers are the bar a Phase-2 model has to clear. The literature
canonical claim is "Burgio RF beats DeepCpf1 by 14-50%" — we want to know
our DeepCpf1's actual Spearman on Kim 2018 endogenous data so we know
what 14-50% means in absolute terms.

Run `tools/fetch_cas12a_data.py` first to populate `data/cas12a/kim_2018/`.

Usage:
    python tools/benchmark_cas12a_baselines.py
"""
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from crisprware.scorers import deepcpf1, enpam_gb

KIM_DIR = os.path.join(REPO_ROOT, "data", "cas12a", "kim_2018")

# (csv filename, label, indel-frequency column name)
EVAL_SHEETS = [
    ("hek_plasmid.csv",  "HEK-plasmid (55 endo)",  "indel_freq_subtracted_pct"),
    ("hct_plasmid.csv",  "HCT-plasmid (66 endo)",  "indel_freq_subtracted_pct"),
    ("hek_lenti.csv",    "HEK-lenti endogenous",   "indel_freq_endogenous_pct"),
    ("hek_lenti.csv",    "HEK-lenti synthetic",    "indel_freq_synthetic_pct"),
    ("kleinstiver_2016.csv", "Kleinstiver 2016 (22 endo)",   "indel_freq_pct"),
    ("chari_2017.csv",       "Chari 2017 (18 endo)",         "indel_freq_pct"),
    ("kim_2016.csv",         "Kim 2016 (10 endo)",           "indel_freq_pct"),
    # Held-out HT sets for reference (synthetic, not endogenous)
    ("ht_2.csv",         "HT 2 (2963 synthetic, held-out)",  "indel_freq_subtracted_pct"),
    ("ht_3.csv",         "HT 3 (1251 synthetic, held-out)",  "indel_freq_subtracted_pct"),
]


def _filter_valid(df, indel_col):
    df = df.copy()
    df["context"] = df["context"].astype(str)
    mask = df["context"].str.len() == 34
    df = df[mask].reset_index(drop=True)
    df[indel_col] = pd.to_numeric(df[indel_col], errors="coerce")
    df = df.dropna(subset=[indel_col]).reset_index(drop=True)
    return df


def _score(seqs, scorer_module, model):
    return np.array(scorer_module.predict(list(seqs), model=model), dtype=float)


def _metrics(scores, truth):
    mask = np.isfinite(scores) & np.isfinite(truth)
    if mask.sum() < 3:
        return float("nan"), float("nan"), int(mask.sum())
    rho, _ = spearmanr(scores[mask], truth[mask])
    r, _ = pearsonr(scores[mask], truth[mask])
    return float(rho), float(r), int(mask.sum())


def main():
    print("Loading models...")
    deep_model = deepcpf1.load_model()
    enpam_model = enpam_gb.load_model()

    print(f"\n{'Sheet':40s}  {'N':>5s}   {'Spearman':>10s}  {'Pearson':>9s}")
    print(f"{'':40s}  {'':>5s}   {'DC | eP':>10s}  {'DC | eP':>9s}")
    print("-" * 80)

    rows = []
    for csv_name, label, indel_col in EVAL_SHEETS:
        path = os.path.join(KIM_DIR, csv_name)
        if not os.path.exists(path):
            print(f"  SKIP {label}: missing {path}")
            continue

        df = pd.read_csv(path)
        df = _filter_valid(df, indel_col)
        if df.empty:
            print(f"  SKIP {label}: no valid rows after filter")
            continue

        truth = df[indel_col].to_numpy(dtype=float)
        dc_scores = _score(df["context"], deepcpf1, deep_model)
        ep_scores = _score(df["context"], enpam_gb, enpam_model)
        dc_rho, dc_r, n_dc = _metrics(dc_scores, truth)
        ep_rho, ep_r, n_ep = _metrics(ep_scores, truth)

        print(f"  {label:38s}  {len(df):>5d}   "
              f"{dc_rho:>4.2f} | {ep_rho:>4.2f}   "
              f"{dc_r:>4.2f} | {ep_r:>4.2f}")
        rows.append({
            "sheet": label,
            "n": len(df),
            "deepcpf1_spearman": dc_rho,
            "enpam_gb_spearman": ep_rho,
            "deepcpf1_pearson": dc_r,
            "enpam_gb_pearson": ep_r,
        })

    summary = pd.DataFrame(rows)
    out_path = os.path.join(REPO_ROOT, "data", "cas12a", "baseline_metrics.csv")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    summary.to_csv(out_path, index=False)
    print(f"\nWrote summary to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
