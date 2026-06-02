#!/usr/bin/env python3
"""Fetch and parse Kim et al. 2018 Supplementary Table 1 into per-sheet CSVs.

The XLSX is hosted on Springer Nature's static CDN under the standard
"Springer Nature" supplementary licensing (free download, derivative-reuse
ambiguous). We download it once into `data/cas12a/raw/` and parse the
relevant sheets into clean tidy CSVs in `data/cas12a/kim_2018/`.

The `data/` tree is gitignored — every checkout regenerates locally.

Usage:
    python tools/fetch_cas12a_data.py
        [--force]               # re-download even if cached
        [--out-dir <path>]      # default: <repo>/data/cas12a

Source: Kim et al. 2018, Nat Biotechnol 36:239,
    https://www.nature.com/articles/nbt.4061
    Supplementary Table 1:
    https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.4061/MediaObjects/41587_2018_BFnbt4061_MOESM39_ESM.xlsx
"""

import argparse
import os
import sys
import urllib.request

ST1_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fnbt.4061/MediaObjects/"
    "41587_2018_BFnbt4061_MOESM39_ESM.xlsx"
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (sheet name in XLSX, output CSV filename, kind)
#   kind = "ht"   -> high-throughput synthetic; has 50-bp, 34-bp, 20-bp guide,
#                    background/Cpf1/subtracted indel frequencies and read counts.
#   kind = "endo" -> endogenous loci; cell line + chromosomal position + 34-bp
#                    context + indel frequency + chromatin accessibility.
#   kind = "endo_lenti" -> same as "endo" but ALSO has a synthetic-site
#                    indel frequency column (HEK-lenti only).
#   kind = "ref"  -> third-party reference set (Kleinstiver/Chari/Kim 2016);
#                    same shape as "endo" with Source column.
SHEETS = [
    ("Data set HT 1-1", "ht_1-1.csv", "ht"),
    ("Data set HT 1-2", "ht_1-2.csv", "ht"),
    ("Data set HT 2", "ht_2.csv", "ht"),
    ("Data set HT 3", "ht_3.csv", "ht"),
    ("Data set HEK-lenti", "hek_lenti.csv", "endo_lenti"),
    ("Data set HEK-plasmid", "hek_plasmid.csv", "endo"),
    ("Data set HCT-plasmid", "hct_plasmid.csv", "endo"),
    ("Kleinstiver 2016", "kleinstiver_2016.csv", "ref"),
    ("Chari 2017", "chari_2017.csv", "ref"),
    ("Kim 2016", "kim_2016.csv", "ref"),
]


def fetch(url, dest, force=False):
    if os.path.exists(dest) and not force:
        print(f"  cached: {dest} ({os.path.getsize(dest):,} bytes)")
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, dest)
    print(f"  -> {dest} ({os.path.getsize(dest):,} bytes)")


def parse_sheet(xlsx_path, sheet_name, kind):
    import pandas as pd

    # Row 0 is the title, row 1 has the per-column descriptions, data begins row 2.
    # We use header=1 to lift the descriptions, then rename to tidy column names.
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name, header=1)
    df = df.dropna(how="all").reset_index(drop=True)

    if kind == "ht":
        rename = {
            df.columns[0]: "context_50bp",
            df.columns[1]: "context",
            df.columns[2]: "guide_20bp",
            df.columns[3]: "indel_freq_background_pct",
            df.columns[4]: "indel_reads_background",
            df.columns[5]: "total_reads_background",
            df.columns[6]: "indel_freq_cpf1_pct",
            df.columns[7]: "indel_reads_cpf1",
            df.columns[8]: "total_reads_cpf1",
            df.columns[9]: "indel_freq_subtracted_pct",
        }
    elif kind == "endo":
        rename = {
            df.columns[0]: "cell_line",
            df.columns[1]: "target_number",
            df.columns[2]: "chromosomal_position_hg19",
            df.columns[3]: "context",
            df.columns[4]: "indel_freq_subtracted_pct",
            df.columns[5]: "chromatin_accessible",
        }
    elif kind == "endo_lenti":
        rename = {
            df.columns[0]: "cell_line",
            df.columns[1]: "target_number",
            df.columns[2]: "chromosomal_position_hg19",
            df.columns[3]: "context",
            df.columns[4]: "indel_freq_endogenous_pct",
            df.columns[5]: "indel_freq_synthetic_pct",
            df.columns[6]: "chromatin_accessible",
        }
    elif kind == "ref":
        rename = {
            df.columns[0]: "source",
            df.columns[1]: "cell_line",
            df.columns[2]: "chromosomal_position_hg19",
            df.columns[3]: "context",
            df.columns[4]: "indel_freq_pct",
            df.columns[5]: "chromatin_accessible",
        }
    else:
        raise ValueError(f"unknown kind {kind!r}")

    df = df.rename(columns=rename)
    df = df[list(rename.values())]
    # Strip any whitespace from sequence columns
    for c in ("context", "context_50bp", "guide_20bp"):
        if c in df.columns:
            df[c] = df[c].astype(str).str.strip()
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--force", action="store_true")
    p.add_argument("--out-dir", default=os.path.join(REPO_ROOT, "data", "cas12a"))
    args = p.parse_args()

    raw_dir = os.path.join(args.out_dir, "raw")
    kim_dir = os.path.join(args.out_dir, "kim_2018")
    os.makedirs(kim_dir, exist_ok=True)

    xlsx = os.path.join(raw_dir, "kim_2018_st1.xlsx")
    fetch(ST1_URL, xlsx, force=args.force)

    print("\nParsing sheets:")
    summary = []
    for sheet_name, out_csv, kind in SHEETS:
        df = parse_sheet(xlsx, sheet_name, kind)
        dest = os.path.join(kim_dir, out_csv)
        df.to_csv(dest, index=False)
        ctx_lens = df["context"].dropna().str.len().value_counts().to_dict() if "context" in df.columns else {}
        summary.append((out_csv, len(df), ctx_lens))
        print(f"  {out_csv:24s} rows={len(df):>6d}  context_len={ctx_lens}")

    print(f"\nWrote {len(summary)} CSVs to {kim_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
