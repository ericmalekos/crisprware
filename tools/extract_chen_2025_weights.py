#!/usr/bin/env python3
"""One-time extraction of Chen 2025 enseq-DeepCpf1 (and optionally the 23
seq-DeepCpf1variants) PyTorch state_dicts into framework-agnostic .npz files.

The runtime crisprware.scorers.enseq_deepcpf1 module loads these .npz files
directly with numpy, so the main env never has to depend on PyTorch.

Source: Chen et al. 2025, Nat Commun 16:3022 (DOI 10.1038/s41467-025-57150-9),
Code Ocean capsule https://codeocean.com/capsule/9398276/tree/v1, source code
under MIT (Yankang Wu, 2024), data under CC0-1.0.

Usage:
    git clone https://git.codeocean.com/capsule-9398276.git
    python tools/extract_chen_2025_weights.py
        [--capsule <path>]       # default: ./capsule-9398276
        [--variants]             # also extract the 23 per-variant models
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEIGHTS_DIR = os.path.join(REPO_ROOT, "crisprware", "scorers", "weights", "chen_2025")


# Maps the EnCas12a (sequence-only) state_dict keys to friendly numpy names.
# Layer order matches the model code (model.py:EnCas12a).
KEY_MAP_ENCAS12A = {
    "layer1.0.weight":    "conv1_w",   # (128, 4, 5) PyTorch -> (5, 4, 128) TF
    "layer1.0.bias":      "conv1_b",   # (128,)
    "conv1ds.0.0.weight": "conv2_w",   # (128, 128, 5) -> (5, 128, 128)
    "conv1ds.0.0.bias":   "conv2_b",   # (128,)
    "layer2.0.weight":    "fc1_w",     # (128, 3968) -> (3968, 128)
    "layer2.0.bias":      "fc1_b",     # (128,)
    "fcs.0.0.weight":     "fc2_w",     # (128, 128) -> (128, 128)
    "fcs.0.0.bias":       "fc2_b",     # (128,)
    "out.0.weight":       "out_w",     # (1, 128) -> (128, 1)
    "out.0.bias":         "out_b",     # (1,)
}


def convert(pth_path, npz_path):
    """Load a PyTorch state_dict, transpose to TF layout, save as .npz."""
    import numpy as np
    import torch

    sd = torch.load(pth_path, map_location="cpu", weights_only=True)
    out = {}
    for k_pt, k_tf in KEY_MAP_ENCAS12A.items():
        if k_pt not in sd:
            raise KeyError(f"{pth_path}: missing key {k_pt}")
        arr = sd[k_pt].numpy()
        if k_tf.endswith("_w") and arr.ndim == 3:
            # Conv1D kernel: PyTorch (out, in, k) -> TF (k, in, out)
            arr = arr.transpose(2, 1, 0)
        elif k_tf.endswith("_w") and arr.ndim == 2:
            # Linear weight: PyTorch (out, in) -> TF (in, out)
            arr = arr.transpose(1, 0)
        out[k_tf] = arr.astype(np.float32)

    os.makedirs(os.path.dirname(npz_path), exist_ok=True)
    np.savez(npz_path, **out)
    print(f"  {os.path.basename(pth_path):40s} -> {os.path.relpath(npz_path, REPO_ROOT)}  "
          f"({sum(a.nbytes for a in out.values()):,} bytes)")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--capsule", default=os.path.join(REPO_ROOT, "capsule-9398276"),
                   help="Path to the cloned Code Ocean capsule")
    p.add_argument("--variants", action="store_true",
                   help="Also extract the 23 per-variant seq-DeepCpf1variants models")
    args = p.parse_args()

    pth_root = os.path.join(args.capsule, "code", "EnDeepCpf1", "trained_model")
    if not os.path.isdir(pth_root):
        print(f"ERROR: {pth_root} not found. Clone the capsule first:", file=sys.stderr)
        print("  git clone https://git.codeocean.com/capsule-9398276.git", file=sys.stderr)
        return 1

    # Main enseq-DeepCpf1
    main_pth = os.path.join(pth_root, "EnDeepCpf1_trained_model.pth")
    convert(main_pth, os.path.join(WEIGHTS_DIR, "enseq_deepcpf1.npz"))

    if args.variants:
        variants_dir = os.path.join(pth_root, "Cas12a_variants")
        out_dir = os.path.join(WEIGHTS_DIR, "variants")
        for fn in sorted(os.listdir(variants_dir)):
            if not fn.endswith(".trained_model.pth"):
                continue
            variant = fn.replace(".trained_model.pth", "").replace(" ", "_")
            convert(
                os.path.join(variants_dir, fn),
                os.path.join(out_dir, f"{variant}.npz"),
            )

    print(f"\nWrote weights to {WEIGHTS_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
