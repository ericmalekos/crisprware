#!/usr/bin/env python3
"""Compare two builds of the same AsCas12a track on CONTENT, not byte layout.

Byte equality is the wrong test. The merge reassigns every `_offset` in chunk order,
so the same guide's off-target list sits at a different byte in each details file, and
the bigBed's R-tree and BGZF block boundaries follow the write order. Two correct
builds of one genome therefore differ in size.

Columns are graded, because they do not all carry the same guarantee:

* HARD: coordinates, sequence, PAM, uniqueness, specificity, `score`, and every
  guide's resolved off-target list. These are exact functions of the genome and the
  CFD matrix. Any difference is a defect.
* SOFT: the cleavage scores and everything downstream of them (their percentiles, the
  efficiency ramp colours, the hover). Keras inference is not batch-invariant, and
  `compute_deepcpf1_scores` slices the guide list into 100,000-sequence calls whose
  phase depends on where the list starts, so a chunked run batches a given guide
  differently from an unchunked one. A score sitting within a hair of a rounding
  boundary can land on either side. Allowed, bounded, and reported.

Joins on (chrom, chromStart, strand, guideSeq), never row-by-row: `sort -k1,1 -k2,2n`
is not a total order for guides sharing a chromStart, and a positional comparison
reports misalignment as if it were a content difference.

Usage:  compare_tracks.py <track_a> <track_b> <enzyme> [soft_ppm_limit]
"""
import os
import re
import subprocess
import sys
import tempfile

HARD = ("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart",
        "thickEnd", "guideSeq", "pam", "unique")
# specificity columns are matched by prefix, since the name is enzyme-specific
HARD_PREFIX = ("AsCas12a_specificity", "enCas12a_spec", "_specColor")
SOFT_SUFFIX = ("Perc",)
SOFT_EXACT = ("itemRgb", "_crisprDeepCpf1Color", "_mouseOver")


def _as_fields(path):
    return [m.group(1) for line in open(path)
            for m in [re.match(r"\s+(?:string|uint|char|bigint)\S*\s+(\w+);", line)] if m]


def _to_bed(bb, out):
    subprocess.run(["bigBedToBed", bb, out], check=True)


def _load(path, key_idx):
    d = {}
    with open(path, "rb") as fh:
        for row in fh:
            f = row.rstrip(b"\n").split(b"\t")
            d[b"|".join(f[i] for i in key_idx)] = f
    return d


def _details_md5(bed, det, key_idx, off_idx):
    """md5 over `key<TAB>off-target line`, sorted. Resolves `_offset` on each side so
    the comparison is of content rather than of byte position."""
    pairs = []
    with open(bed, "rb") as fh:
        for row in fh:
            f = row.rstrip(b"\n").split(b"\t")
            o = int(f[off_idx])
            if o:
                pairs.append((o, b"|".join(f[i] for i in key_idx)))
    pairs.sort(key=lambda t: t[0])
    tmp = tempfile.NamedTemporaryFile(delete=False, dir=os.path.dirname(bed), suffix=".det")
    i = pos = 0
    p = subprocess.Popen(["bgzip", "-dc", "-@", "4", det], stdout=subprocess.PIPE,
                         bufsize=1024 * 1024)
    with tmp:
        for line in p.stdout:
            while i < len(pairs) and pairs[i][0] == pos:
                tmp.write(pairs[i][1] + b"\t" + line)
                i += 1
            pos += len(line)
    p.stdout.close()
    if p.wait() != 0:
        raise RuntimeError(f"bgzip -dc failed for {det}")
    if i != len(pairs):
        raise RuntimeError(f"{i} of {len(pairs)} offsets resolved in {det}")
    s = subprocess.run(f"LC_ALL=C sort -S4G -T {os.path.dirname(bed)} {tmp.name} | md5sum",
                       shell=True, capture_output=True, text=True, check=True)
    os.unlink(tmp.name)
    return s.stdout.split()[0], len(pairs)


def main(a_dir, b_dir, enzyme, ppm_limit="100"):
    ppm_limit = float(ppm_limit)
    bb, det, asf = f"{enzyme}.bb", f"{enzyme}_crisprDetails.tab.gz", f"{enzyme}Targets.as"
    hard_fail = soft_only = 0

    a_as, b_as = os.path.join(a_dir, asf), os.path.join(b_dir, asf)
    ma = subprocess.run(["md5sum", a_as], capture_output=True, text=True).stdout.split()[0]
    mb = subprocess.run(["md5sum", b_as], capture_output=True, text=True).stdout.split()[0]
    print(f"1. autoSql   A {ma}\n             B {mb}\n   {'PASS identical' if ma == mb else 'FAIL differs'}")
    if ma != mb:
        hard_fail += 1

    fields = _as_fields(a_as)
    key_idx = [fields.index(f) for f in ("chrom", "chromStart", "strand", "guideSeq")]
    off_idx = fields.index("_offset")

    work = os.path.join(b_dir, "cmp_work")
    os.makedirs(work, exist_ok=True)
    ab, bb_ = os.path.join(work, "a.bed"), os.path.join(work, "b.bed")
    print("2. bigBed -> BED")
    _to_bed(os.path.join(a_dir, bb), ab)
    _to_bed(os.path.join(b_dir, bb), bb_)

    A, B = _load(ab, key_idx), _load(bb_, key_idx)
    print(f"   A {len(A)} guides, B {len(B)} guides")
    only_a, only_b = set(A) - set(B), set(B) - set(A)
    if only_a or only_b:
        print(f"   FAIL {len(only_a)} guides only in A, {len(only_b)} only in B")
        hard_fail += 1
    else:
        print("   PASS same guide set")

    print("3. per-column comparison (keyed join)")
    per_field = {}
    examples = {}
    for k, fa in A.items():
        fb = B.get(k)
        if fb is None:
            continue
        for i in range(len(fa)):
            if i == off_idx or fa[i] == fb[i]:
                continue
            per_field[i] = per_field.get(i, 0) + 1
            examples.setdefault(i, (k.decode(), fa[i].decode(), fb[i].decode()))
    n = max(len(A), 1)
    if not per_field:
        print("   PASS every column identical on every guide")
    for i, cnt in sorted(per_field.items(), key=lambda t: -t[1]):
        name = fields[i]
        is_hard = (name in HARD or any(name.startswith(p) for p in HARD_PREFIX))
        is_soft = (name in SOFT_EXACT or any(name.endswith(s) for s in SOFT_SUFFIX)
                   or name in ("deepcpf1", "ascas12aDeepcpf1", "enpamGb", "encas12aHf1"))
        tag = "HARD" if is_hard else ("soft" if is_soft else "UNCLASSIFIED")
        k, va, vb = examples[i]
        ppm = 1e6 * cnt / n
        print(f"   {tag:13} {name:24} {cnt:8} rows ({ppm:8.3f} ppm)  {k}  {va!r} -> {vb!r}")
        if is_hard:
            hard_fail += 1
        elif is_soft:
            soft_only = max(soft_only, ppm)
        else:
            hard_fail += 1  # an unclassified column is treated as hard, deliberately

    print("4. every guide's resolved off-target list")
    da, na = _details_md5(ab, os.path.join(a_dir, det), key_idx, off_idx)
    db, nb = _details_md5(bb_, os.path.join(b_dir, det), key_idx, off_idx)
    print(f"   A {da}  ({na} records)\n   B {db}  ({nb} records)")
    if da == db:
        print("   PASS identical off-target lists")
    else:
        print("   FAIL off-target lists differ")
        hard_fail += 1

    for p in (ab, bb_):
        os.path.exists(p) and os.unlink(p)
    print()
    if hard_fail:
        print(f"VERDICT: FAIL  {hard_fail} hard difference(s)")
        return 1
    if soft_only > ppm_limit:
        print(f"VERDICT: FAIL  cleavage drift {soft_only:.3f} ppm exceeds the {ppm_limit} ppm limit")
        return 1
    print(f"VERDICT: PASS  no hard differences; worst cleavage drift {soft_only:.3f} ppm "
          f"(limit {ppm_limit})")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
