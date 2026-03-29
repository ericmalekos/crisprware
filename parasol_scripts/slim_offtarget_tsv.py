#!/usr/bin/env python3
"""
Slim down a scored off-target TSV by streaming out only the columns
needed for downstream merging (combine_on_off_tsvs.py).

Drops the large per-locus columns (offTargets_loci, offTargets_loci_seq)
and otCount, which can dominate file size and memory. Output is gzipped
by default.

Usable as:
  - A standalone CLI:
        python slim_offtarget_tsv.py -i big.reduced -o slim.tsv.gz
        python slim_offtarget_tsv.py -i big.reduced -o slim.tsv --no-gzip

  - A library function from other scripts:
        from parasol_scripts.slim_offtarget_tsv import slim_offtarget_tsv
        slim_offtarget_tsv("big.reduced", "slim.tsv.gz")
"""
from __future__ import annotations

import argparse
import gzip
import sys

# Columns that combine_on_off_tsvs.py actually needs from the off-target file.
DEFAULT_KEEP_COLS = [
    "contig",
    "target",
    "0-mismatch",
    "1-mismatch",
    "2-mismatch",
    "3-mismatch",
    "4-mismatch",
    "TTTN_enCas12a",
    "TTTV_enCas12a",
    "TTTN_enCas12a_aggregated_score",
    "TTTV_enCas12a_aggregated_score",
    "TTTN_2xNLS-Cas12a",
    "TTTV_2xNLS-Cas12a",
    "TTTN_2xNLS-Cas12a_aggregated_score",
    "TTTV_2xNLS-Cas12a_aggregated_score",
    "unique-TTTV",
    "unique-TTTN",
]

# Mismatch columns whose entries contain sequences that can be stripped.
# Format: "SEQUENCE_chrom:pos:strand" -> "chrom:pos:strand"
_MISMATCH_COLS = {"0-mismatch", "1-mismatch", "2-mismatch", "3-mismatch", "4-mismatch"}


def _strip_sequences_from_mismatch(field: str) -> str:
    """Strip leading sequences from mismatch bucket entries.

    'TTTAGTGAAG_chrI:57761:-,TTTGATGATG_chrX:293146:+'
    becomes 'chrI:57761:-,chrX:293146:+'
    """
    if not field or field == "":
        return field
    parts = []
    for entry in field.split(","):
        entry = entry.strip()
        if not entry:
            continue
        if "_" in entry:
            parts.append(entry.rsplit("_", 1)[1])
        else:
            parts.append(entry)
    return ",".join(parts)


def slim_offtarget_tsv(
    input_path: str,
    output_path: str,
    keep_cols: list[str] | None = None,
    gzip_output: bool = True,
) -> str:
    """
    Stream *input_path* line-by-line, writing only *keep_cols* to
    *output_path*.  Only one source line is in memory at a time, so
    the large dropped columns never accumulate.

    Parameters
    ----------
    input_path : str
        Scored off-target TSV (output of score_flashfry_cfd.py).
    output_path : str
        Destination path.  If *gzip_output* is True and the path does
        not already end with '.gz', '.gz' is appended.
    keep_cols : list[str] | None
        Columns to retain.  Columns not present in the file are
        silently skipped.  Defaults to DEFAULT_KEEP_COLS.
    gzip_output : bool
        Gzip the output (default True).

    Returns
    -------
    str
        The actual path written (may have .gz appended).
    """
    if keep_cols is None:
        keep_cols = DEFAULT_KEEP_COLS

    if gzip_output and not output_path.endswith(".gz"):
        output_path += ".gz"

    # -- Peek at header to resolve column indices --
    with open(input_path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")

    header_set = set(header)
    keep_present = [c for c in keep_cols if c in header_set]
    if not keep_present:
        raise ValueError(
            f"None of the requested columns found in {input_path}. "
            f"File columns: {header}"
        )

    col_indices = [header.index(c) for c in keep_present]
    slim_header = "\t".join(keep_present) + "\n"

    # -- Identify which output positions are mismatch columns (need sequence stripping) --
    mismatch_positions = {j for j, name in enumerate(keep_present) if name in _MISMATCH_COLS}

    # -- Stream rows, writing only the slim columns --
    opener = gzip.open(output_path, "wt") if gzip_output else open(output_path, "w")
    with opener as dst, open(input_path, "r") as src:
        dst.write(slim_header)
        next(src)  # skip source header
        for line in src:
            fields = line.rstrip("\n").split("\t")
            out_fields = []
            for j, i in enumerate(col_indices):
                val = fields[i] if i < len(fields) else ""
                if j in mismatch_positions:
                    val = _strip_sequences_from_mismatch(val)
                out_fields.append(val)
            dst.write("\t".join(out_fields) + "\n")

    return output_path


def main():
    ap = argparse.ArgumentParser(
        description="Slim a scored off-target TSV to only the columns needed for merging."
    )
    ap.add_argument(
        "-i", "--input", required=True,
        help="Input scored off-target TSV (from score_flashfry_cfd.py).",
    )
    ap.add_argument(
        "-o", "--output", required=True,
        help="Output path. '.gz' is appended automatically unless --no-gzip.",
    )
    ap.add_argument(
        "--no-gzip", action="store_true", default=False,
        help="Write plain text instead of gzip-compressed output.",
    )
    ap.add_argument(
        "--keep-cols", nargs="+", default=None,
        help="Override which columns to keep (space-separated). "
             "Defaults to the standard set needed by combine_on_off_tsvs.py.",
    )
    args = ap.parse_args()

    try:
        out = slim_offtarget_tsv(
            input_path=args.input,
            output_path=args.output,
            keep_cols=args.keep_cols,
            gzip_output=not args.no_gzip,
        )
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
