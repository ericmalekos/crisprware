#!/usr/bin/env python3
"""One-time conversion: re-save the enPAM+GB joblib so it loads across all
sklearn versions from 0.21 through current (tested 1.0, 1.5+).

Two cross-version incompatibilities the original 0.21.2 pickle hits on
modern sklearn:

  1. `sklearn.ensemble._gb_losses.LeastSquaresError` (referenced by
     `GradientBoostingRegressor.loss_`) was REMOVED in sklearn 1.2.
  2. `sklearn.tree._tree.Tree.__setstate__` in sklearn 1.3+ rejects the
     7-field node dtype, requiring an 8th `missing_go_to_left` (u1) field
     added in 1.3.

This script fixes both by re-saving the joblib with:

  - `loss_` swapped to `crisprware.scorers.enpam_gb.PortableSquaredLoss`,
    a self-contained stub. Eliminates the `_gb_losses` reference.
  - Each `_tree.Tree` re-pickled via a custom reducer pointing at
    `crisprware.scorers.enpam_gb._make_tree`, which tries the original
    `__setstate__` first and pads node dtype on ValueError. The pickle
    bytes still carry the original node dtype, so old sklearn loads
    fine; new sklearn falls through to the pad-and-retry path.

Net: one joblib, every sklearn. The legacy module-rename shim in
`enpam_gb._install_sklearn_legacy_aliases()` is still needed because the
top-level GBR class is referenced as `sklearn.ensemble.gradient_boosting.
GradientBoostingRegressor` -- shim aliases that to the modern module.

Usage:
    python tools/portable_enpam_gb_weights.py

Must be run in an env where the *current* joblib loads -- i.e. sklearn
<= 1.1 (so `_gb_losses` exists), or any sklearn if you've already
applied this conversion once and are just refreshing.
"""

import copyreg
import os
import sys

import numpy as np
import joblib

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

from crisprware.scorers import enpam_gb  # noqa: E402


def _reduce_tree(tree):
    """Pickle-time reducer pointing trees at `_make_tree` for adaptive load."""
    cls, args, state = tree.__reduce__()
    return enpam_gb._make_tree, (args, state)


def main() -> int:
    weights_path = os.path.join(REPO_ROOT, "crisprware", "scorers", "weights", "enpam_gb.joblib")
    if not os.path.exists(weights_path):
        print(f"ERROR: {weights_path} not found", file=sys.stderr)
        return 1

    print(f"Loading {weights_path} ...")
    model = enpam_gb.load_model(weights_path)
    print(f"  loss_  : {type(model.loss_).__module__}.{type(model.loss_).__name__}")

    t0 = model.estimators_[0, 0].tree_
    nodes_dtype = t0.__getstate__()["nodes"].dtype
    print(f"  tree nodes dtype fields: {nodes_dtype.names}")
    print(f"  n estimators           : {model.estimators_.shape}")

    seqs = [
        "ACATTTTACTTTTTCAAAATTGTTTTCATGCTAA",
        "AGGTTTTACAACCGCCCAGTGCGTCTACGTCACA",
        "TTTTTTTAGTGAAGCTTCTAGATATTTGGCGGGT",
    ]
    feat = enpam_gb._featurize_array(seqs)
    before = model.predict(feat)

    if not isinstance(model.loss_, enpam_gb.PortableSquaredLoss):
        model.loss_ = enpam_gb.PortableSquaredLoss()
        print(f"  loss_ -> {type(model.loss_).__name__}")

    from sklearn.tree._tree import Tree

    copyreg.pickle(Tree, _reduce_tree)
    print("  installed copyreg reducer for sklearn.tree._tree.Tree")

    after = model.predict(feat)
    diff = float(np.max(np.abs(before - after)))
    print(f"  max abs diff in-mem (pre-dump): {diff:.2e}")
    if diff > 0:
        print("ERROR: predictions diverged in memory", file=sys.stderr)
        return 1

    joblib.dump(model, weights_path)
    print(f"\nWrote portable joblib to {weights_path} ({os.path.getsize(weights_path):,} bytes)")

    with open(weights_path, "rb") as f:
        blob = f.read()
    refs_gb_losses = b"_gb_losses" in blob
    refs_LSE = b"LeastSquaresError" in blob
    refs_make_tree = b"_make_tree" in blob
    print("\nVerification (new pickle bytes):")
    print(f"  references sklearn.ensemble._gb_losses : {refs_gb_losses}")
    print(f"  references LeastSquaresError           : {refs_LSE}")
    print(f"  references crisprware._make_tree       : {refs_make_tree}")
    if refs_gb_losses or refs_LSE or not refs_make_tree:
        print("WARNING: pickle bytes failed expectation", file=sys.stderr)
        return 1

    print("\nRound-trip check ...")
    reload_model = enpam_gb.load_model(weights_path)
    rt = reload_model.predict(feat)
    rt_diff = float(np.max(np.abs(before - rt)))
    print(f"  max abs diff round-trip: {rt_diff:.2e}")
    if rt_diff > 1e-6:
        print("ERROR: round-trip predictions diverged", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
