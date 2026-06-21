#!/usr/bin/env python3
"""Merge per-chromosome/per-chunk Cas12a UCSC tracks (score_guides --ucscgb) into one.

Thin CLI wrapper around ``crisprware.ucsc_track.merge_ucscgb_tracks`` -- the SAME merge
the native ``score_guides --ucscgb_chunk_max`` path uses, so the standalone fan-out and
the built-in option produce identical output. Each input track dir has:

  guides.bed            bed9+ ; LAST field = _offset = uncompressed byte offset into
                        that piece's crisprDetails.tab (0 = non-unique sentinel).
  crisprDetails.tab.gz  bgzipped; line 0 = header, each later line = one guide's
                        "<mismatchCounts>\\t<offtargetList>".
  cas12aTargets.as      autoSql (identical across pieces).

The off-target details carry no guide key -- the BED's _offset is the only link. The merge
streams each piece's details, maps every data line's old byte-offset -> its new offset in
the merged file, rewrites each piece's BED _offset, then sort + bedToBigBed and bgzip +
index. One shared header sits at byte 0 (the 0 sentinel). Merge order is arbitrary.
"""

import argparse
import os
import sys

# allow running as a standalone script (scripts/ucsc_gb_scripts/ -> repo root)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from crisprware import ucsc_track  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracks", nargs="+", required=True, help="per-piece track dirs (order arbitrary)")
    ap.add_argument("--out", required=True, help="output merged track dir")
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--as-file", required=True, help="cas12aTargets.as (any piece's copy)")
    ap.add_argument("--bedtobigbed", default="bedToBigBed")
    ap.add_argument("--sort-mem", default="8G")
    ap.add_argument("--sort-threads", default="8")
    ap.add_argument("--keep-bed", action="store_true", help="keep guides.sorted.bed (for validation/debug)")
    args = ap.parse_args()

    ucsc_track.merge_ucscgb_tracks(
        track_dirs=args.tracks,
        out=args.out,
        chrom_sizes=args.chrom_sizes,
        as_file=args.as_file,
        bedtobigbed=args.bedtobigbed,
        sort_mem=args.sort_mem,
        sort_threads=args.sort_threads,
        keep_bed=args.keep_bed,
    )


if __name__ == "__main__":
    main()
