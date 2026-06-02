#!/usr/bin/env python3
"""Generate 1000 random NGG-flanked 20-mer guides from a reference genome.

Emits two files in the current working directory:
  - random_1000.fasta — one >id / 23-bp record per guide. Consumed by
    FlashFry's `discover` and crispr-ots's `enumerate --queries`.
  - random_1000.kmers.csv — one row per guide in GuideScan2's
    `--kmers-file` shape. Consumed by GuideScan2 and crispr-ots's
    `enumerate --kmers-file`. Single column whose header is the literal
    string `id,sequence,pam,chromosome,position,sense`.

Usage: make_guides.py <reference.fa>

Seed 42 + sorted-by-discovery-order sampling makes the output
reproducible regardless of Python version.
"""
import gzip
import random
import re
import sys
from pathlib import Path

if len(sys.argv) != 2:
    sys.stderr.write("usage: make_guides.py <reference.fa>\n")
    sys.exit(1)

GENOME = Path(sys.argv[1])
N_GUIDES = 1000
SEED = 42

random.seed(SEED)


def iter_contigs(path):
    """Yield (name, sequence) pairs from a FASTA file. .gz auto-detected."""
    name = None
    chunks = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip().upper())
    if name is not None:
        yield name, "".join(chunks)


NGG_RE = re.compile(r"(?=([ACGT]{20}[ACGT]GG))")


def find_ngg_sites(seq):
    """Return list of (start, 23-bp string) for every NGG-flanked 20-mer."""
    return [(m.start(), m.group(1)) for m in NGG_RE.finditer(seq)]


def main():
    candidates = []
    for name, seq in iter_contigs(GENOME):
        sys.stderr.write(f"  scanning {name} ({len(seq):,} bp)... ")
        sites = find_ngg_sites(seq)
        sys.stderr.write(f"{len(sites):,} NGG sites\n")
        for start, s23 in sites:
            candidates.append((name, start, s23))

    sys.stderr.write(f"\nTotal NGG sites: {len(candidates):,}\n")

    # Sample N_GUIDES, sorted by discovery order so output is stable.
    chosen_idx = sorted(random.sample(range(len(candidates)), N_GUIDES))

    fa = Path("random_1000.fasta").open("w")
    csv = Path("random_1000.kmers.csv").open("w")
    csv.write("id,sequence,pam,chromosome,position,sense\n")
    for idx in chosen_idx:
        chrom, start, s23 = candidates[idx]
        protospacer = s23[:20]
        pam = s23[20:]
        guide_id = f"{chrom}:{start}:+"
        fa.write(f">{guide_id}\n{s23}\n")
        # GuideScan2 wants 1-indexed position (its --kmers-file convention).
        csv.write(f"{guide_id},{protospacer},NGG,{chrom},{start + 1},+\n")
    fa.close()
    csv.close()
    sys.stderr.write(
        f"Wrote {N_GUIDES} guides to random_1000.fasta and random_1000.kmers.csv\n"
    )


if __name__ == "__main__":
    main()
