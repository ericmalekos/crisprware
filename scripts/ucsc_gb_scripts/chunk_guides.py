#!/usr/bin/env python3
"""Split per-chromosome guide BEDs into <=--chunk-max-guide chunks for parallel
`score_guides --ucscgb` scoring.

Input: a directory of per-chromosome BEDs (one ``<chrom>.bed`` each, a header line
followed by data rows) -- e.g. produced by splitting a genome-wide `generate_guides`
BED on column 1. Each output chunk keeps the header line. Writes
``<out-dir>/chunks/<chrom>_p<NN>.bed`` and ``<out-dir>/chunks/manifest.txt`` (one chunk
name per line), suitable for a throttled SLURM job array:

    python chunk_guides.py --bed-dir bed_by_chrom --out-dir . --chunk-max 3000000 \\
        --exclude chr21 chr22
    N=$(wc -l < chunks/manifest.txt)
    sbatch --array=0-$((N-1))%8 score_array.sbatch
    # in the array task: NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" chunks/manifest.txt)

Chunking is purely positional, so the per-chunk tracks recombine with merge_tracks.py
in any order (it reassigns crisprDetails.tab offsets and re-sorts the bigBed).
"""
import argparse
import glob
import math
import os


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed-dir", required=True, help="directory of per-chromosome <chrom>.bed files (header + rows)")
    ap.add_argument("--out-dir", required=True, help="output dir; chunks/ + manifest.txt are written under it")
    ap.add_argument("--chunk-max", type=int, default=3_000_000, help="max guides per chunk [default: 3,000,000]")
    ap.add_argument(
        "--include", nargs="*", default=None, help="explicit chrom names (default: every <chrom>.bed in --bed-dir)"
    )
    ap.add_argument("--exclude", nargs="*", default=[], help="chrom names to skip (e.g. already-scored chromosomes)")
    args = ap.parse_args()

    if args.include:
        chroms = list(args.include)
    else:
        chroms = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(args.bed_dir, "*.bed")))
    skip = set(args.exclude)
    chroms = [c for c in chroms if c not in skip]
    if not chroms:
        raise SystemExit("no chromosomes to chunk (check --bed-dir / --include / --exclude)")

    chunk_dir = os.path.join(args.out_dir, "chunks")
    os.makedirs(chunk_dir, exist_ok=True)
    manifest = []
    for chrom in chroms:
        bed = os.path.join(args.bed_dir, f"{chrom}.bed")
        with open(bed) as fh:
            header = fh.readline()
            n = sum(1 for _ in fh)
        nchunks = max(1, math.ceil(n / args.chunk_max))
        size = math.ceil(n / nchunks)
        with open(bed) as fh:
            fh.readline()  # skip header
            part, written, out = 0, 0, None
            for line in fh:
                if out is None or written == size:
                    if out:
                        out.close()
                    name = f"{chrom}_p{part:02d}"
                    out = open(os.path.join(chunk_dir, f"{name}.bed"), "w")
                    out.write(header)
                    manifest.append(name)
                    part += 1
                    written = 0
                out.write(line)
                written += 1
            if out:
                out.close()
        print(f"{chrom}: {n} guides -> {nchunks} chunks (~{size} each)", flush=True)

    manifest_path = os.path.join(chunk_dir, "manifest.txt")
    with open(manifest_path, "w") as mf:
        mf.write("\n".join(manifest) + "\n")
    print(f"TOTAL {len(manifest)} chunks -> {manifest_path}")


if __name__ == "__main__":
    main()
