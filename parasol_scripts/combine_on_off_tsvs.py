#!/usr/bin/env python3
from __future__ import annotations
import argparse
import os
import sys
import tempfile
import pandas as pd
import numpy as np

def _fmt_pct_raw(raw_s: pd.Series, pct_s: pd.Series) -> pd.Series:
    """
    Return strings like: '92% (0.6782)' — percentile rounded to nearest integer.
    If percentile is missing, return the raw value as-is.
    """
    out = []
    for raw, pct in zip(raw_s, pct_s):
        try:
            f = float(pct)
            out.append(f"{int(round(f))}% ({raw})")
        except Exception:
            out.append(str(raw))
    return pd.Series(out, index=raw_s.index, dtype="object")

def _fmt_mouseover(row) -> str:
    """
    'EnCas12A Spec: #, DeepCPF1: #, EnCas12a: #'
    where # are percentile integers (no % sign).
    """
    def gi(col):
        try:
            return str(int(round(float(row.get(col, "nan")))))
        except Exception:
            return ""
    return (
        f"EnCas12A Spec: {gi('TTTV_enCas12a_aggregated_score_pct')}, "
        f"DeepCpf1: {gi('deepcpf1_score_pct')}, "
        f"EnCas12a: {gi('enpam_gb_score_pct')}"
    )


def _has_any_entries(cell: str) -> bool:
    if pd.isna(cell):
        return False
    toks = [t.strip() for t in str(cell).split(",")]
    return any(t for t in toks)

def _to_float(val):
    try:
        return float(val)
    except Exception:
        return float("nan")

def choose_itemRgb(row) -> str:
    """
    Hierarchy:
      1) If 0-mismatch has any entries -> "150,150,150"
      2) Else if 1- or 2-mismatch has any entries -> "120,120,120"
      3) Else if TTTV_enCas12a_aggregated_score_pct <= 55 -> "81,81,80"
      4) Else color by enpam_gb_score_pct:
           0–50%   -> "230,106,113"
           50–75%  -> "251,243,34"
           75–100% -> "105,181,20"
    """
    # 1)
    if _has_any_entries(row.get("0-mismatch", "")):
        return "150,150,150"
    # 2)
    if _has_any_entries(row.get("1-mismatch", "")) or _has_any_entries(row.get("2-mismatch", "")):
        return "120,120,120"
    # 3)
    tttn_pct = _to_float(row.get("TTTV_enCas12a_aggregated_score_pct", "nan"))
    if not pd.isna(tttn_pct) and tttn_pct <= 55:
        return "81,81,80"
    # 4)
    dcpf_pct = _to_float(row.get("enpam_gb_score_pct", "nan"))
    if pd.isna(dcpf_pct):
        # sensible fallback if percentile missing
        return "150,150,150"
    if dcpf_pct <= 50:
        return "230,106,113"
    elif dcpf_pct <= 75:
        return "251,243,34"
    else:
        return "105,181,20"


def _to_numeric_series(s: pd.Series) -> pd.Series:
    """Coerce to float; treat 'Null'/'', etc. as NaN."""
    return pd.to_numeric(s.replace({"Null": np.nan, "": np.nan}), errors="coerce")

def add_percentile_columns(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Add <col>_pct columns (0–100). Higher value => higher percentile."""
    for c in cols:
        if c in df.columns:
            x = _to_numeric_series(df[c])
            df[c + "_pct"] = (x.rank(pct=True, method="average") * 100).round(2)
    return df

def plot_hist_ecdf(series: pd.Series, title: str, out_png: str):
    """Save a histogram with ECDF overlay (PNG). Skips if series is all NaN."""
    x = _to_numeric_series(series).dropna()
    if x.empty:
        return
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6,4))
    # histogram
    ax.hist(x.values, bins=50, density=True, alpha=0.6)
    ax.set_xlabel(title)
    ax.set_ylabel("Density")
    # ECDF (overlay on twin axis to keep scales readable)
    ax2 = ax.twinx()
    xs = np.sort(x.values)
    ys = np.arange(1, xs.size + 1) / xs.size
    ax2.plot(xs, ys, linewidth=2)
    ax2.set_ylabel("ECDF")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

def load_tsv(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t", dtype=str)
    except Exception as e:
        print(f"ERROR: failed to read {path}: {e}", file=sys.stderr)
        sys.exit(1)

def pick_series(df: pd.DataFrame, name: str, fallback: str = "") -> pd.Series:
    """Return df[name] if present, else a Series of fallback strings of appropriate length."""
    if name in df.columns:
        return df[name]
    return pd.Series([fallback] * len(df), index=df.index, dtype="object")


# --------------------------------------------------------------------------
# Pass-1 helpers: parse one on-target file, merge with one off-target file
# --------------------------------------------------------------------------

REQUIRED_ON = {"chr", "start", "stop", "strand", "enpam_gb_score", "deepcpf1_score"}
REQUIRED_OFF = {"contig", "unique-TTTV", "unique-TTTN",
                "0-mismatch", "1-mismatch", "2-mismatch", "3-mismatch", "4-mismatch"}


def parse_on_target(on_df: pd.DataFrame) -> pd.DataFrame:
    """Parse the comma-packed 4th column and add contig/guideSeq/pam columns."""
    if on_df.shape[1] < 4:
        print("ERROR: on-target TSV must have at least 4 columns.", file=sys.stderr)
        sys.exit(1)

    meta = on_df.iloc[:, 3].astype(str)
    parts = meta.str.split(",", n=5, expand=True)
    if parts.shape[1] != 6:
        print("ERROR: 4th column of on-target TSV must split into 6 items: "
              "id,sequence,pam,chromosome,position,sense", file=sys.stderr)
        sys.exit(1)

    on_df["contig"]   = parts[0].str.strip()
    on_df["guideSeq"] = parts[1].str.strip()
    on_df["pam"]      = parts[2].str.strip()
    return on_df


def validate_columns(df: pd.DataFrame, required: set[str], label: str):
    missing = required.difference(df.columns)
    if missing:
        print(f"ERROR: {label} TSV missing required columns: "
              f"{', '.join(sorted(missing))}", file=sys.stderr)
        sys.exit(1)


def validate_contigs(on_df: pd.DataFrame, off_df: pd.DataFrame,
                     on_path: str, off_path: str):
    """Verify that on-target and off-target contig sets match exactly."""
    on_contigs = set(on_df["contig"].dropna().unique())
    off_contigs = set(off_df["contig"].dropna().unique())

    if on_contigs == off_contigs:
        return

    only_on = on_contigs - off_contigs
    only_off = off_contigs - on_contigs

    msgs = []
    if only_on:
        examples = sorted(only_on)[:5]
        msgs.append(f"  {len(only_on)} contigs only in on-target: {', '.join(examples)}{'...' if len(only_on) > 5 else ''}")
    if only_off:
        examples = sorted(only_off)[:5]
        msgs.append(f"  {len(only_off)} contigs only in off-target: {', '.join(examples)}{'...' if len(only_off) > 5 else ''}")

    raise ValueError(
        f"Contig mismatch between\n"
        f"  on-target:  {on_path}\n"
        f"  off-target: {off_path}\n" +
        "\n".join(msgs)
    )


def merge_pair(on_df: pd.DataFrame, off_df: pd.DataFrame) -> pd.DataFrame:
    """Inner-join one on-target and one off-target DataFrame on contig."""
    return on_df.merge(off_df, on="contig", how="inner", suffixes=("_on", "_off"))


# --------------------------------------------------------------------------
# Completed-log helpers
# --------------------------------------------------------------------------

def load_completed_pairs(log_path: str | None) -> set[tuple[str, str]]:
    """Load set of (on_path, off_path) pairs already completed."""
    if not log_path or not os.path.exists(log_path):
        return set()
    completed = set()
    with open(log_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                completed.add((parts[0], parts[1]))
    return completed


def append_completed(log_path: str, on_path: str, off_path: str):
    """Append a successfully merged pair to the completed log."""
    with open(log_path, "a") as fh:
        fh.write(f"{on_path}\t{off_path}\n")


# --------------------------------------------------------------------------
# Pass-2: read the intermediate, compute percentiles, build BED + mismatch
# --------------------------------------------------------------------------

def build_outputs(merged: pd.DataFrame, args):
    """From the complete merged DataFrame, compute percentiles and write outputs."""
    pam_from_off = merged["target"].astype(str).str.slice(0, 4)

    score_cols = [
        "enpam_gb_score",
        "deepcpf1_score",
        "TTTN_enCas12a_aggregated_score",
        "TTTV_enCas12a_aggregated_score",
    ]
    merged = add_percentile_columns(merged, score_cols)

    if args.plot_dists:
        base, _ = os.path.splitext(args.bed_out)
        for c in score_cols:
            if c in merged.columns:
                out_png = f"{base}.{c}.png"
                plot_hist_ecdf(merged[c], c, out_png)
                print(f"Saved {out_png}")

    itemRgb_series = merged.apply(choose_itemRgb, axis=1)

    # --- Map specificity columns
    TTTV_enCas12a_spec = pick_series(merged, "TTTV_enCas12a_aggregated_score")
    TTTV_enCas12a_spec = np.where(TTTV_enCas12a_spec == "", pick_series(merged, "TTTV_enCas12a"), TTTV_enCas12a_spec)
    TTTN_enCas12a_spec = pick_series(merged, "TTTN_enCas12a_aggregated_score")
    TTTN_enCas12a_spec = np.where(TTTN_enCas12a_spec == "", pick_series(merged, "TTTN_enCas12a"), TTTN_enCas12a_spec)

    _disp_deepcpf1 = _fmt_pct_raw(merged["deepcpf1_score"], merged["deepcpf1_score_pct"])
    _disp_enpam    = _fmt_pct_raw(merged["enpam_gb_score"],  merged["enpam_gb_score_pct"])

    _disp_TTTV_enc = _fmt_pct_raw(
        pd.Series(TTTV_enCas12a_spec, index=merged.index),
        merged.get("TTTV_enCas12a_aggregated_score_pct")
    )
    _disp_TTTN_enc = _fmt_pct_raw(
        pd.Series(TTTN_enCas12a_spec, index=merged.index),
        merged.get("TTTN_enCas12a_aggregated_score_pct")
    )

    _mouseOver = merged.apply(_fmt_mouseover, axis=1)

    # --- Adjust start/stop by strand
    start_raw = pd.to_numeric(merged["start"], errors="coerce").fillna(0).astype(int)
    stop_raw  = pd.to_numeric(merged["stop"],  errors="coerce").fillna(0).astype(int)
    strand    = merged["strand"].astype(str)

    start_adj = np.where(strand == "+", start_raw - 4, start_raw)
    start_adj = np.maximum(start_adj, 0)
    stop_adj  = np.where(strand == "-", stop_raw + 4, stop_raw)

    thick_start = np.where(strand == "+", start_adj + 4, start_adj)
    thick_end   = np.where(strand == "+", stop_adj,      stop_adj - 4)

    # --- Mismatch TSV (compact: counts + coordinates) ---
    # Written first so we can capture byte offsets for the BED file.
    mm_cols = ["0-mismatch", "1-mismatch", "2-mismatch", "3-mismatch", "4-mismatch"]
    have_cols = [c for c in mm_cols if c in merged.columns]

    counts_list = []
    coords_list = []
    for _, row in merged[have_cols].iterrows():
        bucket_counts = []
        all_coords = []
        for col in have_cols:
            cell = row.get(col, "")
            if pd.isna(cell) or not str(cell).strip():
                bucket_counts.append(0)
                continue
            entries = [e.strip() for e in str(cell).split(",") if e.strip()]
            bucket_counts.append(len(entries))
            for entry in entries:
                # Formats supported:
                #   SEQ_chrom:pos:strand[:score]  (full, from .reduced)
                #   chrom:pos:strand[:score]      (slim, sequences stripped)
                if "_" in entry:
                    coord = entry.rsplit("_", 1)[1]
                else:
                    coord = entry
                parts = coord.split(":")
                if len(parts) >= 3:
                    chrom, pos, strand_part = parts[0], parts[1], parts[2]
                    score = parts[3] if len(parts) >= 4 else "0"
                    all_coords.append(f"{chrom};{pos}{strand_part};{score}")
        counts_list.append(",".join(str(c) for c in bucket_counts))
        coords_list.append("|".join(all_coords))

    offsets = []
    with open(args.mismatch_out, "w") as fh:
        fh.write("_mismatchCounts\t_crisprOfftargets\n")
        for counts_str, coords_str in zip(counts_list, coords_list):
            offsets.append(fh.tell())  # byte offset of the start of this line
            fh.write(f"{counts_str}\t{coords_str}\n")

    # --- BED9+ output ---
    bed_df = pd.DataFrame({
        "chr":            merged["chr"],
        "start":          start_adj.astype(int),
        "stop":           stop_adj.astype(int),
        "name":           " ",
        "score":          0,
        "strand":         strand,
        "thickStart":     thick_start.astype(int),
        "thickEnd":       thick_end.astype(int),
        "itemRgb":        itemRgb_series,
        "guideSeq":       merged["guideSeq"],
        "pam":            pam_from_off,
        "deepcpf1_score": _disp_deepcpf1,
        "enpam_gb_score": _disp_enpam,
        "TTTV_enCas12a_specificity": _disp_TTTV_enc,
        "TTTN_enCas12a_specificity": _disp_TTTN_enc,
        "unique-TTTV":    merged["unique-TTTV"],
        "unique-TTTN":    merged["unique-TTTN"],
        "_mouseOver":     _mouseOver,
        "_offset":        offsets,
    })

    bed_df = bed_df.rename(columns={"chr": "#chr"})
    bed_df.to_csv(args.bed_out, sep="\t", index=False)

    if args.minbed:
        from os.path import splitext
        root, ext = splitext(args.bed_out)
        min_path = f"{root}.bed9{ext or '.tsv'}"
        bed9_df = bed_df.iloc[:, :9]
        bed9_df.to_csv(min_path, sep="\t", index=False)
        print(f"Wrote minimal BED9 to {min_path}")


def main():
    ap = argparse.ArgumentParser(
        description="Merge off-target and on-target TSVs; output BED9+ and a mismatch TSV."
    )
    inp = ap.add_mutually_exclusive_group(required=True)
    inp.add_argument("--file-pairs", "-f",
                     help="Two-column TSV (no header) with on-target path in column 1 "
                          "and off-target path in column 2, one pair per line.")
    inp.add_argument("--offtargets", "-x",
                     help="Off-target TSV(s) from score_flashfry_cfd.py "
                          "(comma-separated for multiple files, matched 1:1 with --ontargets)")
    ap.add_argument("--ontargets", "-y",
                    help="On-target TSV(s) (comma-separated for multiple files, "
                         "matched 1:1 with --offtargets). Required when using -x.")
    ap.add_argument("--bed-out", "-b", required=True, help="Output BED9+ TSV path")
    ap.add_argument("--mismatch-out", "-m", required=True, help="Output mismatch TSV path")
    ap.add_argument("--minbed", action="store_true",
                help="Also write a minimal BED (first 9 columns) to a separate file.")
    ap.add_argument("--plot-dists", action="store_true",
                help="If set, save a histogram + CDF PNG for each score column.")
    ap.add_argument("--completed-log", "-c",
                    help="Path to a TSV log that tracks successfully merged pairs. "
                         "On resume, pairs already in this log are skipped. "
                         "Created automatically if it does not exist.")

    args = ap.parse_args()

    # --- Resolve file pairs from either --file-pairs or -x/-y ---
    if args.file_pairs:
        on_paths = []
        off_paths = []
        with open(args.file_pairs) as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) != 2:
                    print(f"ERROR: {args.file_pairs} line {lineno}: expected 2 tab-separated "
                          f"columns, got {len(parts)}", file=sys.stderr)
                    sys.exit(1)
                on_paths.append(parts[0].strip())
                off_paths.append(parts[1].strip())
    else:
        if not args.ontargets:
            print("ERROR: --ontargets/-y is required when using --offtargets/-x.",
                  file=sys.stderr)
            sys.exit(1)
        off_paths = [p.strip() for p in args.offtargets.split(",") if p.strip()]
        on_paths  = [p.strip() for p in args.ontargets.split(",") if p.strip()]

    if len(off_paths) != len(on_paths):
        print(f"ERROR: number of off-target files ({len(off_paths)}) does not match "
              f"number of on-target files ({len(on_paths)}). "
              f"Files must be provided in matched pairs.", file=sys.stderr)
        sys.exit(1)

    if not off_paths:
        print("ERROR: no input files provided.", file=sys.stderr)
        sys.exit(1)

    # --- Load completed pairs for resume support ---
    completed = load_completed_pairs(args.completed_log)
    if completed:
        print(f"Loaded {len(completed)} previously completed pairs from {args.completed_log}")

    # ------------------------------------------------------------------
    # Single pair: fast path — no intermediate file needed
    # ------------------------------------------------------------------
    if len(off_paths) == 1:
        on_path, off_path = on_paths[0], off_paths[0]

        if (on_path, off_path) in completed:
            print(f"Pair already completed, skipping: {off_path} + {on_path}")
        else:
            off_df = load_tsv(off_path)
            on_df  = parse_on_target(load_tsv(on_path))

            validate_columns(on_df, REQUIRED_ON | {"contig", "guideSeq", "pam"}, "on-target")
            validate_columns(off_df, REQUIRED_OFF, "off-target")
            validate_contigs(on_df, off_df, on_path, off_path)

            n_on, n_off = len(on_df), len(off_df)
            merged = merge_pair(on_df, off_df)
            del on_df, off_df

            print(f"on-target rows: {n_on}")
            print(f"off-target rows: {n_off}")
            print(f"merged rows:    {len(merged)}")

            if len(merged) == 0:
                print("WARNING: merge produced 0 rows. Nothing to write.", file=sys.stderr)

            build_outputs(merged, args)

            if args.completed_log:
                append_completed(args.completed_log, on_path, off_path)
        return

    # ------------------------------------------------------------------
    # Multiple pairs: two-pass approach
    #   Pass 1 — merge each pair sequentially, append to intermediate TSV
    #   Pass 2 — read intermediate (compact), compute percentiles, write BED
    # ------------------------------------------------------------------
    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".tsv", prefix="merged_intermediate_")
    os.close(tmp_fd)

    total_on = 0
    total_off = 0
    total_merged = 0
    header_written = False

    try:
        for i, (off_path, on_path) in enumerate(zip(off_paths, on_paths)):
            if (on_path, off_path) in completed:
                print(f"Pair {i+1}/{len(off_paths)} already completed, skipping: "
                      f"{off_path} + {on_path}")
                continue

            print(f"Processing pair {i+1}/{len(off_paths)}: {off_path} + {on_path}")

            off_df = load_tsv(off_path)
            on_df  = parse_on_target(load_tsv(on_path))

            if not header_written:
                validate_columns(on_df, REQUIRED_ON | {"contig", "guideSeq", "pam"}, "on-target")
                validate_columns(off_df, REQUIRED_OFF, "off-target")

            try:
                validate_contigs(on_df, off_df, on_path, off_path)
            except ValueError as e:
                print(f"\nERROR at pair {i+1}/{len(off_paths)}:\n{e}", file=sys.stderr)
                print(f"\n{total_merged} rows from {i} pairs were merged successfully "
                      f"before this error.", file=sys.stderr)
                if args.completed_log:
                    print(f"Re-run with --completed-log {args.completed_log} to resume "
                          f"from pair {i+1}.", file=sys.stderr)
                sys.exit(1)

            n_on, n_off = len(on_df), len(off_df)
            chunk = merge_pair(on_df, off_df)
            del on_df, off_df

            total_on += n_on
            total_off += n_off
            total_merged += len(chunk)

            print(f"  on={n_on}, off={n_off}, merged={len(chunk)}")

            # Append to intermediate file
            chunk.to_csv(tmp_path, sep="\t", index=False,
                         mode="a", header=not header_written)
            header_written = True
            del chunk

            if args.completed_log:
                append_completed(args.completed_log, on_path, off_path)

        print(f"\nTotal on-target rows: {total_on}")
        print(f"Total off-target rows: {total_off}")
        print(f"Total merged rows:    {total_merged}")

        if total_merged == 0:
            print("WARNING: merge produced 0 rows. Nothing to write.", file=sys.stderr)
            return

        # --- Pass 2: read compact intermediate, compute percentiles, write outputs ---
        merged = pd.read_csv(tmp_path, sep="\t", dtype=str)
        build_outputs(merged, args)

    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


if __name__ == "__main__":
    main()
