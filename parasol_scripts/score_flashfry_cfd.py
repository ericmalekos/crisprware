#!/usr/bin/env python3
from __future__ import annotations
import argparse
import pandas as pd
import re
from typing import List, Dict, Tuple
from Bio import SeqIO
import gzip


OFFTARGET_SEQ_RE = re.compile(r"^([ACGT]+)")  # leading DNA seq
TOKEN_RE = re.compile(r"^(?P<seq>[ACGT]+)_(?P<count>\d+)_(?P<mm>\d+)<(?P<sites>[^>]+)>$")
SITE_RE = re.compile(r"(?P<chrom>[^:<>|,]+):(?P<pos>\d+)\^(?P<strand>[RF])")
CONTIG_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<pos>\d+):(?P<strand>[+-])$")


def _split_zero_mismatch_cell(cell: str) -> list[str]:
    """
    From a '0-mismatch' cell like 'SEQ1_chrA:1:+,SEQ2_chrB:2:-',
    return just the sequences ['SEQ1', 'SEQ2'] (strip everything after first underscore).
    """
    if pd.isna(cell) or not str(cell).strip():
        return []
    seqs = []
    for tok in str(cell).split(","):
        tok = tok.strip()
        if not tok:
            continue
        seq = tok.split("_", 1)[0].strip()  # before the first underscore
        if seq:
            seqs.append(seq)
    return seqs


def _has_ttTV(seqs: list[str]) -> bool:
    """True if ANY sequence starts with TTT followed by A/C/G (i.e., not T)."""
    for s in seqs:
        if len(s) >= 4 and s.startswith("TTT") and s[3] in ("A", "C", "G"):
            return True
    return False


def _has_tttt(seqs: list[str]) -> bool:
    """True if ANY sequence starts with TTTT."""
    return any(s.startswith("TTTT") for s in seqs)


def _expected_canonical_coord_from_contig(contig: str) -> str | None:
    """
    From 'chr:pos:strand' compute the expected 0-mismatch coordinate of the *original* target
    that we want to remove from the 0-mismatch list:
      + : pos-1
      - : pos-21
    Returns 'chr:pos:strand' or None if parsing fails.
    """
    if pd.isna(contig) or not str(contig).strip():
        return None
    m = CONTIG_RE.match(str(contig).strip())
    if not m:
        return None
    chrom = m.group("chrom")
    pos = int(m.group("pos"))
    strand = m.group("strand")
    if strand == "+":
        pos = pos - 1
    else:
        pos = pos - 21
    if pos < 1:
        return None
    return f"{chrom}:{pos}:{strand}"


def _remove_canonical_zero_mismatch(zcell: str, contig: str) -> str:
    """
    Remove the entry whose coord matches the expected canonical coord from contig.
    Supports both formats:
      - New: 'chrom:pos:strand:score,...'
      - Legacy: 'SEQ_chrom:pos:strand,...'
    """
    exp = _expected_canonical_coord_from_contig(contig)
    if not exp or pd.isna(zcell) or not str(zcell).strip():
        return "" if pd.isna(zcell) else str(zcell).strip()

    keep = []
    for tok in str(zcell).split(","):
        tok = tok.strip()
        if not tok:
            continue
        # Extract coord (chrom:pos:strand) from entry
        if "_" in tok:
            # Legacy format: SEQ_chrom:pos:strand
            coord = tok.rsplit("_", 1)[1].strip()
        else:
            # New format: chrom:pos:strand:score — first 3 colon fields
            parts = tok.split(":")
            coord = ":".join(parts[:3]) if len(parts) >= 3 else tok
        if coord != exp:
            keep.append(tok)
    return ",".join(keep)


def count_mismatches(a: str, b: str, trim5: int = 4, cap_len: int = 20) -> int:
    """
    Count mismatches after trimming the first `trim5` bases from BOTH sequences,
    then comparing up to `cap_len` positions (default: 20).
    Safe if sequences aren't exactly 27nt.
    """
    if not a or not b:
        return 999
    a = str(a).strip().upper()[trim5:]
    b = str(b).strip().upper()[trim5:]
    L = min(cap_len, len(a), len(b))
    return sum(1 for i in range(L) if a[i] != b[i])


def bucket_offtargets_by_mismatch(
    target_seq: str,
    loci_seqs_csv: str,
    loci_csv: str,
    lookup: Dict[Tuple[int, str], float] | None = None,
    max_bucket: int = 4,
    cfd_scale: int = 1000,
) -> Tuple[dict[int, list[str]], float, float]:
    """
    Pair offTargets_loci_seq with offTargets_loci (index-aligned), compare to target_seq,
    bucket by mismatch count, and compute per-entry CFD scores.

    Returns:
        (buckets, tttn_sum, tttv_sum)
        buckets: {0: [...], 1: [...], ...} with entries like 'chrI:1031460:+:152'
                 (coord:scaledCFD where scaledCFD = int(CFD * cfd_scale))
        tttn_sum: sum of all raw CFD scores (for aggregated TTTN specificity)
        tttv_sum: sum of raw CFD scores excluding TTTT-prefixed entries (for TTTV)
    """
    buckets = {i: [] for i in range(max_bucket + 1)}
    tttn_sum = 0.0
    tttv_sum = 0.0
    if not loci_seqs_csv or not loci_csv or not target_seq:
        return buckets, tttn_sum, tttv_sum

    seqs = [s for s in map(str.strip, str(loci_seqs_csv).split(",")) if s]
    coords = [c for c in map(str.strip, str(loci_csv).split(",")) if c]

    for s, c in zip(seqs, coords):
        mm = count_mismatches(target_seq, s)
        if 0 <= mm <= max_bucket:
            if lookup is not None:
                cfd = cfd_score_pair(target_seq, s, lookup)
                scaled = int(cfd * cfd_scale)
                buckets[mm].append(f"{c}:{scaled}")
                tttn_sum += cfd
                if not s.startswith("TTTT"):
                    tttv_sum += cfd
            else:
                buckets[mm].append(f"{s}_{c}")
    return buckets, tttn_sum, tttv_sum


def expand_target_from_contig(genome, contig_val: str) -> str:
    """
    Parse 'chr:pos:strand' from the contig column and return a 27-mer:
      '+' strand: [pos+1 .. pos+27]
      '-' strand: [pos-2 .. pos+24] then reverse-complement
    Returns "" if parsing fails, chrom missing, or out-of-bounds.
    """
    if pd.isna(contig_val) or not str(contig_val).strip():
        return ""
    m = CONTIG_RE.match(str(contig_val).strip())
    if not m:
        return ""
    chrom = m.group("chrom")
    pos = int(m.group("pos"))
    strand = m.group("strand")

    if strand == "-":
        pos -= 21
    else:
        pos -= 1
    mer = fetch_27mer(genome, chrom, pos, strand)
    return mer if mer is not None else ""


def _open_fasta(path: str):
    # Text mode for SeqIO.parse

    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def load_genome_biopython(fa_path: str):
    """
    Load FASTA/FASTA.gz using Biopython.
    Returns dict: chrom -> SeqRecord
    """
    with _open_fasta(fa_path) as handle:
        return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))


def fetch_27mer(genome_dict, chrom: str, pos_1based: int, strand: str):
    """
    '+' : 27 nt to the right: [pos+1 .. pos+27]
    '-' : 3 before & 24 after: [pos-2 .. pos+24], then reverse-complement
    Returns uppercase 27-mer string or None if OOB/missing chrom.
    """
    rec = genome_dict.get(chrom)
    if rec is None:
        return None
    n = len(rec)
    if strand == "+":
        start = pos_1based + 1
        end = pos_1based + 27
    else:
        start = pos_1based - 2
        end = pos_1based + 24

    if start < 1 or end > n:
        return None

    window = rec.seq[start - 1 : end]  # Seq object
    if strand == "-":
        window = window.reverse_complement()

    return str(window).upper()


def extract_loci_from_cell(cell: str) -> str:
    """
    Parse an offTargets cell like:
      SEQ_1_4<chrX:16599922^R>,SEQ_2_4<chrIII:6555340^F|chrIII:6543296^F>
    into:
      'chrX:16599922:-,chrIII:6555340:+,chrIII:6543296:+'
    """
    if pd.isna(cell) or not str(cell).strip():
        return ""
    out = []
    for tok in str(cell).split(","):
        tok = tok.strip()
        if not tok:
            continue
        m = TOKEN_RE.match(tok)
        if m:
            sites_field = m.group("sites")
        else:
            # Fallback: just pull anything inside <...>
            angle = re.search(r"<([^>]+)>", tok)
            if not angle:
                continue
            sites_field = angle.group(1)

        for site in sites_field.split("|"):
            site = site.strip()
            m2 = SITE_RE.match(site)
            if not m2:
                continue
            strand = "+" if m2.group("strand") == "F" else "-"
            out.append(f"{m2.group('chrom')}:{m2.group('pos')}:{strand}")
    return ",".join(out)


def parse_offtargets(cell) -> List[str]:
    if pd.isna(cell) or str(cell).strip() == "":
        return []
    seqs = []
    for item in str(cell).split(","):
        m = OFFTARGET_SEQ_RE.match(item.strip())
        if m:
            seqs.append(m.group(1))
    return seqs


def parse_offtargets_with_counts(cell) -> List[Tuple[str, int]]:
    """Parse offTargets cell, returning (sequence, q_i) tuples.

    Each FlashFry token has format: SEQ_COUNT_MM<loci>
    COUNT is q_i — the number of genomic loci where that off-target occurs.
    """
    if pd.isna(cell) or str(cell).strip() == "":
        return []
    results = []
    for item in str(cell).split(","):
        item = item.strip()
        if not item:
            continue
        m = TOKEN_RE.match(item)
        if m:
            results.append((m.group("seq"), int(m.group("count"))))
        else:
            # Fallback: try to extract just the sequence with count=1
            m2 = OFFTARGET_SEQ_RE.match(item)
            if m2:
                results.append((m2.group(1), 1))
    return results


def build_score_lookup(df_score: pd.DataFrame) -> Tuple[str, Dict[Tuple[int, str], float]]:
    required = {"RDA", "Pos", "MM", "avg_percent_active"}
    if not required.issubset(df_score.columns):
        raise ValueError("Score CSV must have columns: RDA, Pos, MM, avg_percent_active")
    matrix_name = str(df_score["RDA"].iloc[0])
    df_score = df_score.copy()
    df_score["Pos"] = df_score["Pos"].astype(int)
    df_score["MM"] = df_score["MM"].astype(str)
    df_score["avg_percent_active"] = df_score["avg_percent_active"].astype(float)
    lookup = {(int(r.Pos), str(r.MM)): float(r.avg_percent_active) for _, r in df_score.iterrows()}
    return matrix_name, lookup


def dna_to_rna_base(b: str) -> str:
    b = b.upper()
    return "U" if b == "T" else b


def cfd_score_pair(target_dna: str, offtarget_dna: str, lookup: Dict[Tuple[int, str], float]) -> float:
    t = target_dna.strip().upper()
    o = offtarget_dna.strip().upper()
    if len(t) < 5 or len(o) < 5:
        return 1.0
    t_sub = t[4:]
    o_sub = o[4:]
    L = min(23, len(t_sub), len(o_sub))
    score = 1.0
    for i in range(L):
        tb, ob = t_sub[i], o_sub[i]
        if tb != ob:
            mm = f"r{dna_to_rna_base(tb)}:d{ob}"
            pos = i + 1
            score *= lookup.get((pos, mm), 1.0)
    return float(score)


def format_scores(scores: List[float]) -> str:
    return ",".join(str(s) for s in scores)


def parse_score_list(cell: str) -> List[float]:
    if pd.isna(cell):
        return []
    vals = []
    for tok in str(cell).split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            vals.append(float(tok))
        except ValueError:
            # e.g., "Null" or non-numeric tokens
            pass
    return vals


def _round_float_str(cell: str, nd: int = 6) -> str:
    """Round a scalar numeric string to nd decimals; keep 'Null' or empty as-is."""
    if cell is None:
        return ""
    s = str(cell).strip()
    if s == "" or s.lower() == "null":
        return s
    try:
        return f"{float(s):.{nd}f}"
    except ValueError:
        return s  # leave non-numeric tokens untouched


def _round_csv_list(cell: str, nd: int = 6) -> str:
    """
    Round a comma-separated list of numeric tokens to nd decimals.
    Keeps 'Null' tokens and empties as-is; trims whitespace.
    """
    if cell is None:
        return ""
    s = str(cell).strip()
    if s == "":
        return s
    parts = []
    for tok in s.split(","):
        tok = tok.strip()
        if tok == "" or tok.lower() == "null":
            parts.append(tok)
            continue
        try:
            parts.append(f"{float(tok):.{nd}f}")
        except ValueError:
            parts.append(tok)  # non-numeric, keep as-is
    return ",".join(parts)


def main():
    ap = argparse.ArgumentParser(description="Append CFD-like scores from mismatch matrices to FlashFry TSV.")
    ap.add_argument(
        "-i", "--input", required=True, help="FlashFry TSV input (needs 'contig', 'target' and 'offTargets')."
    )
    ap.add_argument(
        "-m",
        "--matrices",
        required=True,
        nargs="+",
        help="One or more scoring matrix CSVs (RDA,Pos,MM,avg_percent_active).",
    )
    ap.add_argument("-o", "--output", required=True, help="Output TSV path with appended score columns.")
    ap.add_argument("-g", "--genome", required=True, help="Reference genome in FASTA or FASTA.gz")
    ap.add_argument(
        "--no-gzip",
        action="store_true",
        default=False,
        help="Write the slim output as plain TSV instead of gzip-compressed.",
    )
    args = ap.parse_args()

    # ---- Load & validate input ----
    df = pd.read_csv(args.input, sep="\t", dtype=str)
    required = {"contig", "target", "offTargets", "orientation", "start", "stop"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Input TSV missing required columns: {', '.join(sorted(missing))}")

    # ---- Normalize & filter rows ----
    df["orientation"] = df["orientation"].astype(str).str.strip().str.upper()
    df["start"] = df["start"].astype(str).str.strip()
    df["stop"] = df["stop"].astype(str).str.strip()

    df = df[(df["orientation"] == "FWD") & (df["start"] == "0") & (df["stop"] == "24")].copy()
    df.reset_index(drop=True, inplace=True)

    # After filtering + reset_index
    df = df[(df["orientation"] == "FWD") & (df["start"] == "0") & (df["stop"] == "24")].copy()
    df.reset_index(drop=True, inplace=True)

    # Drop unneeded columns
    df = df.drop(columns=["start", "stop", "context", "overflow", "orientation"], errors="ignore")

    # ---- Parse loci BEFORE any offTargets mutation ----
    df["offTargets_loci"] = df["offTargets"].apply(extract_loci_from_cell)

    # ---- Load genome (Biopython) & extract 27-mers ----
    genome = load_genome_biopython(args.genome)

    # per-locus off-target sequences (comma-separated 27-mers)
    def sequences_from_loci(loci_str: str) -> str:
        if not loci_str:
            return ""
        out = []
        for tok in loci_str.split(","):
            tok = tok.strip()
            if not tok:
                continue
            try:
                chrom, pos, strand = tok.split(":")
                pos = int(pos)
                strand = "+" if strand == "+" else "-"
            except Exception:
                continue
            mer = fetch_27mer(genome, chrom, pos, strand)
            if mer is not None:
                out.append(mer)
        return ",".join(out)

    df["offTargets_loci_seq"] = df["offTargets_loci"].apply(sequences_from_loci)

    # ---- NEW: expand/replace target from contig 27-mer ----
    df["target"] = df["contig"].apply(lambda s: expand_target_from_contig(genome, s))

    # ---- Load scoring matrices ----
    lookups: List[Tuple[str, Dict[Tuple[int, str], float]]] = []
    for matrix_path in args.matrices:
        score_df = pd.read_csv(matrix_path)
        matrix_name, lookup = build_score_lookup(score_df)
        lookups.append((matrix_name, lookup))

    primary_lookup = lookups[0][1] if lookups else None

    # --- Bucket off-targets by mismatch count, scoring each entry with CFD ---
    zero, one, two, three, four = [], [], [], [], []
    agg_sums: Dict[str, Tuple[List[float], List[float]]] = {name: ([], []) for name, _ in lookups}

    for _, row in df.iterrows():
        target_seq = row.get("target", "")
        loci_seqs = row.get("offTargets_loci_seq", "")
        loci_csv = row.get("offTargets_loci", "")

        # Primary bucket + score (first matrix, scores embedded in entries)
        buckets, _, _ = bucket_offtargets_by_mismatch(
            target_seq=target_seq,
            loci_seqs_csv=loci_seqs,
            loci_csv=loci_csv,
            lookup=primary_lookup,
            max_bucket=4,
        )
        zero.append(",".join(buckets[0]))
        one.append(",".join(buckets[1]))
        two.append(",".join(buckets[2]))
        three.append(",".join(buckets[3]))
        four.append(",".join(buckets[4]))

        # Compute aggregation sums for ALL matrices
        for matrix_name, lkp in lookups:
            _, tttn_sum, tttv_sum = bucket_offtargets_by_mismatch(
                target_seq=target_seq,
                loci_seqs_csv=loci_seqs,
                loci_csv=loci_csv,
                lookup=lkp,
                max_bucket=4,
            )
            agg_sums[matrix_name][0].append(tttn_sum)
            agg_sums[matrix_name][1].append(tttv_sum)

    df["0-mismatch"] = zero
    df["1-mismatch"] = one
    df["2-mismatch"] = two
    df["3-mismatch"] = three
    df["4-mismatch"] = four

    df["0-mismatch"] = df.apply(
        lambda r: _remove_canonical_zero_mismatch(r.get("0-mismatch", ""), r.get("contig", "")), axis=1
    )

    # ---- Aggregated scores: 1 / sum(CFD_i) per matrix ----
    for matrix_name, _ in lookups:
        tttn_sums, tttv_sums = agg_sums[matrix_name]
        df[f"TTTN_{matrix_name}_aggregated_score"] = [f"{1.0 / ts:.6f}" if ts > 0 else "" for ts in tttn_sums]
        df[f"TTTV_{matrix_name}_aggregated_score"] = [f"{1.0 / vs:.6f}" if vs > 0 else "" for vs in tttv_sums]

    # ---- Determine if unique in genome ----
    # 0-mismatch entries are now coord:score format. We need the original
    # sequences to determine TTTV/TTTT PAM status, so cross-reference
    # with offTargets_loci_seq and offTargets_loci.
    uniq_ttTV_vals = []
    uniq_tttn_vals = []

    for _, row in df.iterrows():
        cell = row.get("0-mismatch", "")
        if pd.isna(cell) or not str(cell).strip():
            uniq_ttTV_vals.append(True)
            uniq_tttn_vals.append(True)
            continue

        entries = [e.strip() for e in str(cell).split(",") if e.strip()]
        if not entries:
            uniq_ttTV_vals.append(True)
            uniq_tttn_vals.append(True)
            continue

        # Extract coords from 0-mismatch entries (format: chrom:pos:strand:score)
        zero_coords = set()
        for e in entries:
            parts = e.split(":")
            if len(parts) >= 3:
                zero_coords.add(f"{parts[0]}:{parts[1]}:{parts[2]}")

        # Match back to loci sequences to check PAM prefix
        loci_seqs = [s for s in map(str.strip, str(row.get("offTargets_loci_seq", "")).split(",")) if s]
        loci_coords = [c for c in map(str.strip, str(row.get("offTargets_loci", "")).split(",")) if c]

        has_tttv = False
        has_tttt = False
        for s, c in zip(loci_seqs, loci_coords):
            if c in zero_coords:
                if s.startswith("TTTT"):
                    has_tttt = True
                elif len(s) >= 4 and s[:3] == "TTT" and s[3] in "ACG":
                    has_tttv = True

        if has_tttv:
            uniq_ttTV_vals.append(False)
            uniq_tttn_vals.append(False)
        elif has_tttt:
            uniq_ttTV_vals.append(True)
            uniq_tttn_vals.append(False)
        else:
            uniq_ttTV_vals.append(False)
            uniq_tttn_vals.append(False)

    df["unique-TTTV"] = pd.Series(uniq_ttTV_vals, dtype="boolean")
    df["unique-TTTN"] = pd.Series(uniq_tttn_vals, dtype="boolean")

    # ---- Drop unneeded columns before saving ----
    drop_cols = {"offTargets", "offTargets_loci", "offTargets_loci_seq", "otCount"}
    df = df[[c for c in df.columns if c not in drop_cols]]

    # ---- Write output ----
    if args.no_gzip:
        df.to_csv(args.output, sep="\t", index=False)
    else:
        if not args.output.endswith(".gz"):
            args.output += ".gz"
        df.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
