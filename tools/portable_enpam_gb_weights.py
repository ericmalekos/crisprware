#!/usr/bin/env python3
"""One-time conversion: re-save the enPAM+GB joblib with no sklearn-version
deps in its pickle.

The original joblib was pickled with sklearn 0.21.2. Its
`GradientBoostingRegressor.loss_` is an instance of
`sklearn.ensemble._gb_losses.LeastSquaresError`, a class that was removed
entirely in sklearn 1.2. Loading the joblib on a modern sklearn raises
`ModuleNotFoundError: No module named 'sklearn.ensemble._gb_losses'`
even with the module-rename shims in place.

This script:
  1. Loads the existing joblib in the current env (the shim in
     `crisprware.scorers.enpam_gb._install_sklearn_legacy_aliases`
     handles the sklearn 0.22+ module renames -- you need sklearn 1.0.x
     or 1.1.x to run this, since 1.2+ has no `_gb_losses` at all).
  2. Replaces `model.loss_` with a self-contained stub
     (`crisprware.scorers.enpam_gb.PortableSquaredLoss`) that lives in
     our own package, so the re-saved pickle has zero references to
     `sklearn.ensemble._gb_losses`.
  3. Sanity-checks that predictions stay bit-exact (max abs diff = 0).
  4. Re-dumps to `crisprware/scorers/weights/enpam_gb.joblib`.

After this conversion the joblib loads cleanly in any sklearn from 0.22
through current without needing the rename shims either. The shims are
kept for back-compat in case someone is running an environment with a
pre-conversion joblib still on disk.

Usage:
    python tools/portable_enpam_gb_weights.py
"""

import os
import sys

import numpy as np
import joblib

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from crisprware.scorers import enpam_gb


def main() -> int:
    weights_path = os.path.join(
        REPO_ROOT, "crisprware", "scorers", "weights", "enpam_gb.joblib"
    )
    if not os.path.exists(weights_path):
        print(f"ERROR: {weights_path} not found", file=sys.stderr)
        return 1

    print(f"Loading {weights_path} ...")
    model = enpam_gb.load_model(weights_path)
    print(f"  loss_ class: {type(model.loss_).__module__}.{type(model.loss_).__name__}")

    seqs = [
        "ACATTTTACTTTTTCAAAATTGTTTTCATGCTAA",
        "AGGTTTTACAACCGCCCAGTGCGTCTACGTCACA",
        "TTTTTTTAGTGAAGCTTCTAGATATTTGGCGGGT",
    ]
    feat = enpam_gb._featurize_array(seqs)
    before = model.predict(feat)

    model.loss_ = enpam_gb.PortableSquaredLoss()
    print(f"  loss_ class: {type(model.loss_).__module__}.{type(model.loss_).__name__}")

    after = model.predict(feat)
    diff = float(np.max(np.abs(before - after)))
    print(f"  max abs diff vs original predict: {diff:.2e}")
    if diff > 0:
        print("ERROR: predictions diverged after swap", file=sys.stderr)
        return 1

    joblib.dump(model, weights_path)
    print(f"\nWrote portable joblib to {weights_path} ({os.path.getsize(weights_path):,} bytes)")

    with open(weights_path, "rb") as f:
        blob = f.read()
    refs_gb_losses = b"_gb_losses" in blob
    refs_LSE = b"LeastSquaresError" in blob
    print("\nVerification (new pickle bytes):")
    print(f"  references sklearn.ensemble._gb_losses : {refs_gb_losses}")
    print(f"  references LeastSquaresError           : {refs_LSE}")
    if refs_gb_losses or refs_LSE:
        print("WARNING: stale references still in the pickle", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
