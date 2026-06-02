#!/usr/bin/env python3
"""Compare `id -> specificity` across crispr-ots / GuideScan2 / FlashFry.

Usage:
    compare_specs.py <crispr_ots_csv> <gs2_csv> <flashfry_scored> <our_tsv>
                     [kmers_csv]

The first three files are the per-tool outputs. The fourth is our TSV
(FlashFry-convention specificity) used for the FlashFry comparison —
the CSV uses the GS2 convention, so we re-run with --spec-convention
flashfry to get a FlashFry-comparable column. The optional fifth arg
is the kmers CSV used to bridge `id` ↔ `target` (defaults to
`random_1000.kmers.csv` in the current working directory).

Emits a brief diff summary suitable for pasting into BENCHMARKS.md.
"""

import sys


def load_crispr_ots_csv(path):
    """Per-guide id -> specificity. CSV has one row per (guide, hit);
    specificity is repeated, so first occurrence wins."""
    out = {}
    with open(path) as fh:
        header = fh.readline().strip().split(",")
        id_i = header.index("id")
        sp_i = header.index("specificity")
        for line in fh:
            f = line.rstrip("\n").split(",")
            if len(f) <= sp_i:
                continue
            gid = f[id_i]
            if gid not in out:
                out[gid] = float(f[sp_i])
    return out


def load_gs2_csv(path):
    """Same shape as crispr-ots CSV."""
    return load_crispr_ots_csv(path)


def load_flashfry_scored(path):
    """FlashFry scored TSV: key by `target` (23-bp protospacer+PAM),
    value DoenchCFD_specificityscore."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            t_i = header.index("target")
            s_i = header.index("DoenchCFD_specificityscore")
        except ValueError as e:
            sys.stderr.write(f"FlashFry header missing column: {e}\n")
            return {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= s_i:
                continue
            target = f[t_i]
            try:
                out[target] = float(f[s_i])
            except ValueError:
                pass
    return out


def load_our_tsv(path):
    """Our TSV: key by target, value cfd_specificity."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            t_i = header.index("target")
            s_i = header.index("cfd_specificity")
        except ValueError as e:
            sys.stderr.write(f"Our TSV header missing column: {e}\n")
            return {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= s_i:
                continue
            out[f[t_i]] = float(f[s_i])
    return out


def read_id_to_target(kmers_csv):
    """The kmers CSV ties id -> protospacer+pam concat. We use this
    to bridge FlashFry's target-keyed scored output to GS2's id-keyed
    CSV output."""
    out = {}
    with open(kmers_csv) as fh:
        fh.readline()  # header
        for line in fh:
            f = line.strip().split(",")
            if len(f) < 3:
                continue
            gid, proto, pam = f[0], f[1], f[2]
            out[gid] = proto + pam
    return out


def summarize_diff(label, a, b, tol):
    common = set(a) & set(b)
    drift = 0
    worst = 0.0
    for k in common:
        d = abs(a[k] - b[k])
        if d > worst:
            worst = d
        if d > tol:
            drift += 1
    print(f"  {label}: {len(common)} common, max |Δ| = {worst:.3e}, # > {tol:.0e} = {drift}")


def main(argv):
    if len(argv) not in (5, 6):
        sys.stderr.write(__doc__)
        sys.exit(1)
    ours_csv_p, gs2_csv_p, ff_scored_p, ours_tsv_p = argv[1:5]
    kmers_csv_p = argv[5] if len(argv) == 6 else "random_1000.kmers.csv"

    ours_csv = load_crispr_ots_csv(ours_csv_p)
    gs2 = load_gs2_csv(gs2_csv_p)
    ours_tsv = load_our_tsv(ours_tsv_p)
    ff = load_flashfry_scored(ff_scored_p)
    id_to_target = read_id_to_target(kmers_csv_p)

    # Bridge: GS2 keys by id; FlashFry keys by target. Translate the GS2 /
    # crispr-ots CSV results to target-keyed dicts for direct compare.
    def by_target(by_id):
        out = {}
        for gid, sp in by_id.items():
            t = id_to_target.get(gid)
            if t:
                out[t] = sp
        return out

    ours_csv_byT = by_target(ours_csv)
    gs2_byT = by_target(gs2)

    print(f"Loaded {len(ours_csv)} crispr-ots CSV, {len(gs2)} GS2 CSV rows by id")
    print(f"Loaded {len(ff)} FlashFry, {len(ours_tsv)} crispr-ots TSV rows by target")
    print()
    print("Pairwise specificity diff:")
    # GS2 convention: bit-identical (~1e-6 print precision).
    summarize_diff("crispr-ots (GS2 conv) vs GuideScan2", ours_csv_byT, gs2_byT, 1e-5)
    # FlashFry convention: bit-identical to ULP.
    summarize_diff("crispr-ots (FF conv)  vs FlashFry  ", ours_tsv, ff, 1e-9)


if __name__ == "__main__":
    main(sys.argv)
