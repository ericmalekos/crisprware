#!/usr/bin/env python3
"""Phase-2 v1: train a Cas12a on-target activity model on Kim 2018 HT 1-1
and benchmark against the vendored DeepCpf1 / enPAM+GB ports.

Three model variants in parallel:
  - sklearn.ensemble.GradientBoostingRegressor (same class as enPAM+GB,
    apples-to-apples model swap)
  - LightGBM with default-ish hyperparameters
  - XGBoost with default-ish hyperparameters

Features: the 689-dim RS2 set already vendored in
crisprware.scorers.enpam_gb (positional-independent 1-/2-mers,
positional-dependent 1-/2-mers, GC, Tm).

Training:   Kim 2018 HT 1-1 (15,000 AsCas12a synthetic, the set DeepCpf1
            also trained on -- so HT 2 / endogenous remain held out for
            both).
Evaluation: HT 2 (synthetic held-out, 2,963), plus the four endogenous
            sheets and three small third-party reference sets.

Outputs a per-model Spearman/Pearson table side-by-side with the existing
DeepCpf1 / enPAM+GB baseline. The "Phase 2 success bar" is +10 pp over
DeepCpf1's endogenous Spearman.

Run `tools/fetch_cas12a_data.py` first to populate data/cas12a/kim_2018/.

Usage:
    python tools/train_cas12a_v1.py [--quick]
"""

import argparse
import os
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from crisprware.scorers import deepcpf1, enpam_gb

KIM_DIR = os.path.join(REPO_ROOT, "data", "cas12a", "kim_2018")

EVAL_SHEETS = [
    ("hek_plasmid.csv", "HEK-plasmid (endo)", "indel_freq_subtracted_pct"),
    ("hct_plasmid.csv", "HCT-plasmid (endo)", "indel_freq_subtracted_pct"),
    ("hek_lenti.csv", "HEK-lenti endo", "indel_freq_endogenous_pct"),
    ("hek_lenti.csv", "HEK-lenti synthetic", "indel_freq_synthetic_pct"),
    ("kleinstiver_2016.csv", "Kleinstiver 2016", "indel_freq_pct"),
    ("chari_2017.csv", "Chari 2017", "indel_freq_pct"),
    ("kim_2016.csv", "Kim 2016", "indel_freq_pct"),
    ("ht_2.csv", "HT 2 (synthetic)", "indel_freq_subtracted_pct"),
    ("ht_3.csv", "HT 3 (synthetic)", "indel_freq_subtracted_pct"),
]


def _load_dataset(csv_name, indel_col):
    df = pd.read_csv(os.path.join(KIM_DIR, csv_name))
    df["context"] = df["context"].astype(str)
    df = df[df["context"].str.len() == 34].reset_index(drop=True)
    df[indel_col] = pd.to_numeric(df[indel_col], errors="coerce")
    df = df.dropna(subset=[indel_col]).reset_index(drop=True)
    return df


def _spearman(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    return float(spearmanr(a[mask], b[mask])[0])


def _pearson(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    return float(pearsonr(a[mask], b[mask])[0])


def _featurize(seqs):
    return enpam_gb.featurize(seqs).to_numpy()


def train_sklearn_gbr(X, y, quick=False):
    from sklearn.ensemble import GradientBoostingRegressor

    n = 50 if quick else 500
    model = GradientBoostingRegressor(
        n_estimators=n,
        learning_rate=0.05,
        max_depth=3,
        subsample=0.9,
        random_state=42,
        validation_fraction=0.1,
        n_iter_no_change=20,
        tol=1e-4,
    )
    model.fit(X, y)
    return model


def train_lightgbm(X, y, X_val, y_val, quick=False):
    import lightgbm as lgb

    n = 100 if quick else 2000
    model = lgb.LGBMRegressor(
        n_estimators=n,
        learning_rate=0.05,
        num_leaves=31,
        min_data_in_leaf=20,
        feature_fraction=0.9,
        bagging_fraction=0.9,
        bagging_freq=5,
        random_state=42,
        verbose=-1,
    )
    model.fit(
        X,
        y,
        eval_set=[(X_val, y_val)],
        callbacks=[lgb.early_stopping(50, verbose=False)],
    )
    return model


def train_xgboost(X, y, X_val, y_val, quick=False):
    import xgboost as xgb

    n = 100 if quick else 2000
    model = xgb.XGBRegressor(
        n_estimators=n,
        learning_rate=0.05,
        max_depth=6,
        subsample=0.9,
        colsample_bytree=0.9,
        random_state=42,
        early_stopping_rounds=50,
        verbosity=0,
    )
    model.fit(X, y, eval_set=[(X_val, y_val)], verbose=False)
    return model


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--quick", action="store_true", help="Use fewer iterations for fast turnaround")
    args = p.parse_args()

    print("Loading training data (Kim 2018 HT 1-1)...")
    train_df = _load_dataset("ht_1-1.csv", "indel_freq_subtracted_pct")
    print(f"  HT 1-1: {len(train_df)} valid rows")

    print("Loading validation data (Kim 2018 HT 1-2)...")
    val_df = _load_dataset("ht_1-2.csv", "indel_freq_subtracted_pct")
    print(f"  HT 1-2: {len(val_df)} valid rows")

    print("Featurizing (RS2, 689 dims)...")
    t0 = time.time()
    X_train = _featurize(train_df["context"].tolist())
    y_train = train_df["indel_freq_subtracted_pct"].to_numpy(dtype=float)
    print(f"  X_train: {X_train.shape}  ({time.time() - t0:.1f}s)")

    t0 = time.time()
    X_val = _featurize(val_df["context"].tolist())
    y_val = val_df["indel_freq_subtracted_pct"].to_numpy(dtype=float)
    print(f"  X_val: {X_val.shape}  ({time.time() - t0:.1f}s)")

    # --- Train ---
    print("\nTraining sklearn GradientBoostingRegressor...")
    t0 = time.time()
    sk_model = train_sklearn_gbr(X_train, y_train, quick=args.quick)
    print(f"  fit: {time.time() - t0:.1f}s  iters used: {sk_model.n_estimators_}")

    print("Training LightGBM...")
    t0 = time.time()
    lgb_model = train_lightgbm(X_train, y_train, X_val, y_val, quick=args.quick)
    print(f"  fit: {time.time() - t0:.1f}s  iters used: {lgb_model.best_iteration_}")

    print("Training XGBoost...")
    t0 = time.time()
    xgb_model = train_xgboost(X_train, y_train, X_val, y_val, quick=args.quick)
    print(f"  fit: {time.time() - t0:.1f}s  iters used: {xgb_model.best_iteration}")

    # --- Load baselines once ---
    print("\nLoading DeepCpf1 + enPAM+GB baselines...")
    deep_model = deepcpf1.load_model()
    enpam_model = enpam_gb.load_model()

    # --- Evaluate ---
    header = (
        f"\n  {'Sheet':28s}  {'N':>5s}   {'DC':>4s}  {'eP':>4s} | {'sk':>4s}  {'lgb':>4s}  {'xgb':>4s}    (Spearman)"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))
    rows = []
    for csv_name, label, indel_col in EVAL_SHEETS:
        df = _load_dataset(csv_name, indel_col)
        if df.empty:
            continue
        seqs = df["context"].tolist()
        y_true = df[indel_col].to_numpy(dtype=float)
        X_eval = _featurize(seqs)

        dc = np.asarray(deepcpf1.predict(seqs, model=deep_model), dtype=float)
        ep = np.asarray(enpam_gb.predict(seqs, model=enpam_model), dtype=float)
        sk = sk_model.predict(X_eval)
        lg = lgb_model.predict(X_eval)
        xg = xgb_model.predict(X_eval)

        ro = {
            "sheet": label,
            "n": len(df),
            "deepcpf1": _spearman(dc, y_true),
            "enpam_gb": _spearman(ep, y_true),
            "sk_gbr": _spearman(sk, y_true),
            "lightgbm": _spearman(lg, y_true),
            "xgboost": _spearman(xg, y_true),
        }
        rows.append(ro)
        print(
            f"  {label:28s}  {ro['n']:>5d}   "
            f"{ro['deepcpf1']:>4.2f}  {ro['enpam_gb']:>4.2f} | "
            f"{ro['sk_gbr']:>4.2f}  {ro['lightgbm']:>4.2f}  {ro['xgboost']:>4.2f}"
        )

    summary = pd.DataFrame(rows)
    out = os.path.join(REPO_ROOT, "data", "cas12a", "v1_metrics.csv")
    summary.to_csv(out, index=False)
    print(f"\nWrote {out}")

    # --- Phase 2 success-bar table ---
    print("\n  Phase 2 endogenous deltas (model - DeepCpf1 baseline, target >= +0.10):")
    endo = summary[summary["sheet"].str.contains("endo", case=False)]
    for _, r in endo.iterrows():
        for col, name in [("sk_gbr", "sk_gbr"), ("lightgbm", "lgb"), ("xgboost", "xgb")]:
            delta = r[col] - r["deepcpf1"]
            mark = "OK" if delta >= 0.10 else "  "
            print(f"    {r['sheet']:24s}  {name:6s}  {r[col]:>5.2f}  d={delta:+.2f}  {mark}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
