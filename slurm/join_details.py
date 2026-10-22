#!/usr/bin/env python3
"""Emit `key<TAB>details-line` for every guide in a track, so two tracks built by
different paths can be compared on CONTENT rather than on byte layout.

The chunked and unchunked builds cannot have byte-identical details files: the merge
reassigns every `_offset` in chunk order, so the same guide's off-target list sits at a
different byte in each file. `_offset` is therefore the one column that is allowed to
differ, and comparing the files directly would report a difference that means nothing.

Resolving it is a merge join, not a lookup table: the BED's (offset -> key) pairs are
sorted numerically and the decompressed details are streamed once in offset order, so
memory is O(guides sharing one offset) rather than O(guides).

Usage:  join_details.py <guides.bed> <details.tab.gz> <out.tsv>
"""
import subprocess
import sys

KEY_FIELDS = (0, 1, 5, 11)  # chrom, chromStart, strand, guideSeq
OFFSET_FIELD = 20


def main(bed: str, details: str, out: str) -> int:
    # 1. (offset, key) for every guide with a details record. Offset 0 is the
    #    non-unique sentinel pointing at the shared header, which carries no
    #    per-guide content, so those guides have nothing to compare here; the
    #    column-wise check covers them.
    pairs = []
    with open(bed, "rb") as fh:
        for row in fh:
            f = row.rstrip(b"\n").split(b"\t")
            off = int(f[OFFSET_FIELD])
            if off == 0:
                continue
            pairs.append((off, b":".join(f[i] for i in KEY_FIELDS)))
    pairs.sort(key=lambda t: t[0])
    print(f"{bed}: {len(pairs)} guides with a details record", flush=True)

    # 2. stream the details once, in offset order, emitting each needed line.
    n = 0
    i = 0
    with open(out, "wb") as outf:
        p = subprocess.Popen(["bgzip", "-dc", "-@", "4", details],
                             stdout=subprocess.PIPE, bufsize=1024 * 1024)
        pos = 0
        for line in p.stdout:
            while i < len(pairs) and pairs[i][0] == pos:
                outf.write(pairs[i][1] + b"\t" + line)
                i += 1
                n += 1
            pos += len(line)
        p.stdout.close()
        if p.wait() != 0:
            raise RuntimeError(f"bgzip -dc failed for {details}")

    if i != len(pairs):
        unmatched = pairs[i][0] if i < len(pairs) else None
        raise RuntimeError(
            f"{i} of {len(pairs)} offsets resolved; first unresolved offset {unmatched}. "
            "An offset that lands mid-line means the BED and the details disagree."
        )
    print(f"  resolved {n} details lines -> {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
