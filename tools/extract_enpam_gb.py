#!/usr/bin/env python3
"""One-time extraction of the enPAM+GB model from upstream sgrna_modeler.

Run this once in a throwaway env to drop a clean sklearn-only joblib into
`crisprware/scorers/weights/enpam_gb.joblib`. After that, the main crisprware
environment depends only on scikit-learn at runtime — no `sgrna_modeler`,
no joblib download from a third-party URL during normal use.

Usage
-----
    cd <repo root>
    python -m venv /tmp/enpam_extract
    /tmp/enpam_extract/bin/pip install sgrna_modeler scikit-learn joblib
    /tmp/enpam_extract/bin/python tools/extract_enpam_gb.py

Source
------
sgrna_modeler is MIT-licensed (Peter DeWeirdt, Doench lab):
    https://github.com/PeterDeWeirdt/sgrna_modeler
The shipped enPAM+GB model accompanies Luo et al. 2020,
    https://doi.org/10.1038/s41587-020-0600-6

This script must run on Python 3.6+ (the basilisk env that ships with the
crisprScore R package uses Python 3.6).
"""

import os
import sys

try:
    import joblib
    from sgrna_modeler import enzymes as en
    from sgrna_modeler import models as sg
except ImportError as e:
    print(f"ERROR: missing dependency: {e}", file=sys.stderr)
    print(__doc__, file=sys.stderr)
    sys.exit(1)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(REPO_ROOT, "crisprware", "scorers", "weights")
OUT_PATH = os.path.join(OUT_DIR, "enpam_gb.joblib")


def main():
    upstream_path = sg.get_enpam_gb()
    print(f"Upstream model file: {upstream_path}")
    print(f"  size: {os.path.getsize(upstream_path):,} bytes")

    wrapped = sg.SklearnSgrnaModel()
    wrapped.load_model(upstream_path, en.cas12a, "enPAM_GB")
    model = wrapped.model
    print(f"Loaded model: {type(model).__module__}.{type(model).__name__}")
    print(f"  features default (RS2): {wrapped.features}")
    print(
        f"  enzyme: guide_start={en.cas12a['guide_start']} "
        f"guide_length={en.cas12a['guide_length']} "
        f"context_length={en.cas12a['context_length']}"
    )

    if hasattr(model, "n_estimators_"):
        print(f"  n_estimators_: {model.n_estimators_}")
    if hasattr(model, "n_features_in_"):
        print(f"  n_features_in_: {model.n_features_in_}")

    os.makedirs(OUT_DIR, exist_ok=True)
    joblib.dump(model, OUT_PATH)
    print(f"Wrote clean sklearn-only joblib to {OUT_PATH}")
    print(f"  size: {os.path.getsize(OUT_PATH):,} bytes")

    print()
    print("Next: in the main crisprware env, run")
    print("  pytest tests/test_enpam_gb.py -v")
    return 0


if __name__ == "__main__":
    sys.exit(main())
