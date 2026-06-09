#!/usr/bin/env python3
"""Merge per-chromosome Cas12a UCSC tracks (score_guides --ucscgb) into one track.

Each input track dir has:
  guides.bed            bed9+ ; LAST field = _offset = uncompressed byte offset into
                        that chrom's crisprDetails.tab (0 = non-unique sentinel).
  crisprDetails.tab.gz  bgzipped; line 0 = header, each later line = one guide's
                        "<mismatchCounts>\\t<offtargetList>".
  cas12aTargets.as      autoSql (identical across chroms).

The off-target details carry no guide key — the BED's _offset is the only link. So we
stream each chrom's details, map every data line's old byte-offset -> its new offset in
the merged file, rewrite that chrom's BED _offset, then sort + bedToBigBed and
bgzip + index the merged details. One shared header sits at byte 0 (the 0 sentinel).
"""

import argparse
import os
import subprocess
import sys

HEADER = b"_mismatchCounts\t_crisprOfftargets\n"


def zcat_lines(path):
    """Yield raw byte lines (including trailing newline) from a bgzipped file."""
    p = subprocess.Popen(["bgzip", "-dc", "-@", "4", path], stdout=subprocess.PIPE, bufsize=1024 * 1024)
    try:
        for line in p.stdout:
            yield line
    finally:
        p.stdout.close()
        if p.wait() != 0:
            raise RuntimeError(f"bgzip -dc failed for {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracks", nargs="+", required=True, help="per-chrom track dirs, in chrom-sorted order")
    ap.add_argument("--out", required=True, help="output merged track dir")
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--as-file", required=True, help="cas12aTargets.as (any chrom's copy)")
    ap.add_argument("--bedtobigbed", default="bedToBigBed")
    ap.add_argument("--sort-mem", default="8G")
    ap.add_argument("--sort-threads", default="8")
    ap.add_argument("--keep-bed", action="store_true", help="keep guides.sorted.bed (for validation/debug)")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    merged_tab = os.path.join(args.out, "crisprDetails.tab")
    merged_bed = os.path.join(args.out, "guides.merged.bed")

    n_guides = n_unique = 0
    running = 0  # byte offset in the merged uncompressed details
    with open(merged_tab, "wb") as out, open(merged_bed, "wb") as bedout:
        out.write(HEADER)
        running = len(HEADER)
        for tdir in args.tracks:
            det = os.path.join(tdir, "crisprDetails.tab.gz")
            bed = os.path.join(tdir, "guides.bed")
            if not (os.path.exists(det) and os.path.exists(bed)):
                sys.exit(f"missing track files in {tdir}")
            # 1. stream this chrom's details -> old(byte in chrom file) -> new(byte in merged)
            offmap = {}
            chrom_off = 0
            for i, line in enumerate(zcat_lines(det)):
                if i == 0:  # per-chrom header at byte 0 -> not copied (merged has its own at 0)
                    if line != HEADER:
                        sys.exit(f"unexpected details header in {det}: {line!r}")
                    chrom_off += len(line)
                    continue
                offmap[chrom_off] = running
                out.write(line)
                running += len(line)
                chrom_off += len(line)
            # 2. rewrite this chrom's BED with remapped _offset (last column)
            with open(bed, "rb") as bf:
                for row in bf:
                    row = row.rstrip(b"\n")
                    if not row:
                        continue
                    fields = row.split(b"\t")
                    old = int(fields[-1])
                    if old == 0:
                        new = 0
                    else:
                        try:
                            new = offmap[old]
                        except KeyError:
                            sys.exit(f"BED _offset {old} not found in {det} (offset map mismatch)")
                        n_unique += 1
                    fields[-1] = str(new).encode()
                    bedout.write(b"\t".join(fields) + b"\n")
                    n_guides += 1
            offmap.clear()
            print(f"  merged {os.path.basename(tdir)}: cum_guides={n_guides} merged_bytes={running}", flush=True)

    # 3. sort the merged BED and build the bigBed
    sorted_bed = os.path.join(args.out, "guides.sorted.bed")
    subprocess.run(
        f"LC_ALL=C sort -S{args.sort_mem} --parallel={args.sort_threads} -T {args.out} "
        f"-k1,1 -k2,2n {merged_bed} > {sorted_bed}",
        shell=True,
        check=True,
    )
    bb = os.path.join(args.out, "minimumCas12A.bb")
    subprocess.run(
        [args.bedtobigbed, "-tab", "-type=bed9+", f"-as={args.as_file}", sorted_bed, args.chrom_sizes, bb], check=True
    )
    # 4. bgzip + index the merged details (offsets are uncompressed; .gzi handles the seek)
    subprocess.run(["bgzip", "-f", "-@", "8", merged_tab], check=True)  # -> crisprDetails.tab.gz
    subprocess.run(["bgzip", "-r", merged_tab + ".gz"], check=True)  # -> crisprDetails.tab.gz.gzi
    subprocess.run(["cp", args.as_file, os.path.join(args.out, "cas12aTargets.as")], check=True)
    os.remove(merged_bed)
    if not args.keep_bed:
        os.remove(sorted_bed)
    print(f"\nMERGE DONE -> {args.out}\n  guides={n_guides}  unique(with details)={n_unique}  details_bytes={running}")


if __name__ == "__main__":
    main()
