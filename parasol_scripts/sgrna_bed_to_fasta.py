#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import sys
import math


def write_fasta(ids, seqs, output_path):
    with open(output_path, "w") as f:
        for fasta_id, seq in zip(ids, seqs):
            f.write(f">{fasta_id}\n{seq}\n")


def main():
    parser = argparse.ArgumentParser(description="Convert sgRNA BED-style TSV to FASTA.")
    parser.add_argument("input_file", help="Input TSV file with sgRNA annotations.")
    parser.add_argument("output_file", help="Output FASTA file (or prefix if splitting).")
    parser.add_argument("--split", action="store_true", help="Split output by chromosome.")
    parser.add_argument("--flashfry_split", action="store_true", help="Split into FASTAs with ≤1,048,575 entries each.")
    args = parser.parse_args()

    # Read TSV, skipping commented lines
    df = pd.read_csv(args.input_file, sep="\t", comment="#", header=None)

    # Parse ID and sequence
    fasta_ids = df[3].str.split(",", expand=True)[0]
    sequences = df[4].str.slice(start=4, stop=4 + 27)

    df["fasta_id"] = fasta_ids
    df["fasta_seq"] = sequences

    if args.flashfry_split:
        max_per_file = 1000000
        total_chunks = math.ceil(len(df) / max_per_file)
        for i in range(total_chunks):
            chunk_df = df.iloc[i * max_per_file : (i + 1) * max_per_file]
            suffix = f"_{i + 1:06d}.fa"
            out_path = args.output_file + suffix
            write_fasta(chunk_df["fasta_id"], chunk_df["fasta_seq"], out_path)
            print(f"Wrote {len(chunk_df)} entries to {out_path}", file=sys.stderr)

    elif args.split:
        for chrom, group in df.groupby(0):
            chrom_safe = chrom.replace("/", "_")
            out_path = f"{args.output_file}_{chrom_safe}.fa"
            write_fasta(group["fasta_id"], group["fasta_seq"], out_path)
            print(f"Wrote {len(group)} entries to {out_path}", file=sys.stderr)

    else:
        write_fasta(df["fasta_id"], df["fasta_seq"], args.output_file)


if __name__ == "__main__":
    main()
