#!/usr/bin/env python3

import argparse
import os
import pandas as pd


def split_tsv(input_file, n, outdir=None):
    df = pd.read_csv(input_file, sep="\t")
    total_rows = len(df)
    chunk_size = (total_rows + n - 1) // n  # Ceiling division

    # Determine base name and output directory
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        stem = os.path.splitext(os.path.basename(input_file))[0]
        base_prefix = os.path.join(outdir, stem)
    else:
        # Default: keep original behavior (write alongside input file)
        base_prefix = os.path.splitext(input_file)[0]

    for i in range(n):
        start = i * chunk_size
        end = min(start + chunk_size, total_rows)
        chunk_df = df.iloc[start:end]

        if len(chunk_df) == 0:
            continue  # Skip empty chunks

        output_file = f"{base_prefix}.{i + 1}.tsv"
        chunk_df.to_csv(output_file, sep="\t", index=False)

        # Sanity check: remove file if it's empty
        if os.path.exists(output_file) and os.path.getsize(output_file) == 0:
            os.remove(output_file)
        else:
            print(f"Wrote {len(chunk_df)} rows to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Split a TSV file into N parts.")
    parser.add_argument("tsv_file", help="Path to the input TSV file")
    parser.add_argument("-n", "--n", type=int, required=True, help="Number of output files")
    parser.add_argument(
        "-o", "--outdir", default=None, help="Directory to write output files (defaults to input file's directory)"
    )
    args = parser.parse_args()

    split_tsv(args.tsv_file, args.n, args.outdir)


if __name__ == "__main__":
    main()
