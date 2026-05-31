#!/usr/bin/env python3
"""Generate 1000 random TTTN-flanked 20-mer Cas12a guides from a reference.

Cas12a (Cpf1) PAM is at the 5' end: 4-bp `TTTN` precedes the 20-bp
protospacer. Each emitted record is 24 bp:  `TTTN` + protospacer.

Emits two files in cwd:
  - random_1000_cpf1.fasta — `>id` + 24-bp record. Consumed by FlashFry's
    `discover` (Cpf1 mode) and crispr-ots's `enumerate --queries`.
  - random_1000_cpf1.kmers.csv — GuideScan2 / crisprware shape header
    `id,sequence,pam,chromosome,position,sense`. Consumed by crispr-ots
    `enumerate --kmers-file` (GuideScan2 doesn't ship Cas12a support).

Usage: make_guides_cpf1.py <reference.fa[.gz]>
"""
import gzip
import random
import re
import sys
from pathlib import Path

if len(sys.argv) != 2:
    sys.stderr.write("usage: make_guides_cpf1.py <reference.fa[.gz]>\n")
    sys.exit(1)

GENOME = Path(sys.argv[1])
N_GUIDES = 1000
SEED = 42

random.seed(SEED)


def iter_contigs(path):
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


# Lookahead so overlapping matches all fire. 4-bp TTTN PAM + 20-bp protospacer.
TTTN_RE = re.compile(r"(?=(TTT[ACGT][ACGT]{20}))")


def find_tttn_sites(seq):
    return [(m.start(), m.group(1)) for m in TTTN_RE.finditer(seq)]


def main():
    candidates = []
    for name, seq in iter_contigs(GENOME):
        sys.stderr.write(f"  scanning {name} ({len(seq):,} bp)... ")
        sites = find_tttn_sites(seq)
        sys.stderr.write(f"{len(sites):,} TTTN sites\n")
        for start, s24 in sites:
            candidates.append((name, start, s24))

    sys.stderr.write(f"\nTotal TTTN sites: {len(candidates):,}\n")

    chosen_idx = sorted(random.sample(range(len(candidates)), N_GUIDES))

    fa = Path("random_1000_cpf1.fasta").open("w")
    csv = Path("random_1000_cpf1.kmers.csv").open("w")
    csv.write("id,sequence,pam,chromosome,position,sense\n")
    for idx in chosen_idx:
        chrom, start, s24 = candidates[idx]
        pam = s24[:4]                # TTTN
        protospacer = s24[4:]        # 20 bp
        guide_id = f"{chrom}:{start}:+"
        fa.write(f">{guide_id}\n{s24}\n")
        csv.write(f"{guide_id},{protospacer},TTTN,{chrom},{start + 1},+\n")
    fa.close()
    csv.close()
    sys.stderr.write(
        f"Wrote {N_GUIDES} guides to random_1000_cpf1.fasta and random_1000_cpf1.kmers.csv\n"
    )


if __name__ == "__main__":
    main()
