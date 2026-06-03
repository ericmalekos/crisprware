#!/usr/bin/env python3
"""Build a UCSC Genome Browser Cas12a track from crispr-ots scored output.

Consumes the per-guide on-target scores plus the crispr-ots streaming
off-target output (Mode-1 aggregated CSV + Mode-2 per-off-target TSV + the
guide sidecar) and emits the three files pinned by
``cas12a_track_output_spec.md``:

* ``minimumCas12A.bb``        bigBed (bed9+) of on-target guides
* ``crisprDetails.tab.gz``    bgzip per-guide off-target details (+ ``.gzi``)
* ``cas12aTargets.as``        autoSql schema (build input)

Coordinate conventions (validated against hg38, both strands):

* On-target feature spans PAM(4) + spacer(23) = 27 nt. The guide BED's
  ``start``/``stop`` bound the 23-nt spacer on the + strand; the 4-nt PAM is
  added on the 5' side per strand (``+``: left, ``-``: right).
* Off-target ``pos`` (the browser's 24-nt window = 4 PAM + 20 match) maps from
  crispr-ots's 0-based 27-nt-site ``start`` as **``+``: pos = start;
  ``-``: pos = start + 3``** (the 23→20 truncation drops the 3 PAM-distal nt).

The ``_mismatchCounts`` come from Mode-1 (true per-bucket totals, no CFD floor);
the ``_crisprOfftargets`` list comes from Mode-2 (floored at the run's
``--cfd-threshold``), sorted by score and capped per guide for disk control.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

SITE_LEN = 27
SPACER_LEN = 23
PAM_LEN = 4

# On-target display score columns, in track column order. Any subset that is
# actually present in the guide table is emitted.
DISPLAY_SCORE_COLS = ["deepcpf1_score", "enpam_gb_score", "enseq_deepcpf1_score"]

# Track filenames (must not be renamed — the browser config references them).
BB_NAME = "minimumCas12A.bb"
DETAILS_NAME = "crisprDetails.tab"
AS_NAME = "cas12aTargets.as"


# ----------------------------------------------------------------------------
# autoSql
# ----------------------------------------------------------------------------
def build_autosql(score_cols: Sequence[str]) -> str:
    """autoSql for the bigBed: bed9 + guideSeq/pam + on-target/specificity
    display columns + the required ``_mouseOver``/``_offset`` trailer."""
    lines = [
        "table cas12aTargets",
        '"CRISPR Cas12a guides, genome wide (bed9 + extra; off-targets in external bgzip tab file)"',
        "    (",
        '    string  chrom;        "Reference sequence chromosome or scaffold"',
        '    uint    chromStart;   "0-based start of the 27nt PAM+spacer site"',
        '    uint    chromEnd;     "chromStart + 27"',
        '    string  name;         "Guide ID (may be empty)"',
        '    uint    score;        "Score 0-1000 (filter configured 0-100)"',
        '    char[1] strand;       "+ or -"',
        '    uint    thickStart;   "Start of spacer (thick); the 4nt PAM is the thin part"',
        '    uint    thickEnd;     "End of spacer (thick)"',
        '    uint    itemRgb;      "Display color R,G,B"',
        "    string  guideSeq;     \"Spacer 5'->3' (23nt)\"",
        "    string  pam;          \"4nt PAM 5'->3' (e.g. TTTA)\"",
    ]
    for c in score_cols:
        lines.append(f'    string  {c};{" " * max(1, 18 - len(c))}"free-form score (display only)"')
    lines += [
        '    string  TTTV_enCas12a_specificity; "free-form score (display only)"',
        '    string  TTTN_enCas12a_specificity; "free-form score (display only)"',
        '    string  unique_TTTV;               "free-form (display only)"',
        '    string  unique_TTTN;               "free-form (display only)"',
        '    string  _mouseOver;   "Hover label"',
        '    bigint  _offset;      "Uncompressed byte offset into crisprDetails.tab, or 0 if none"',
        "    )",
        "",
    ]
    return "\n".join(lines)


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def _pct(series: pd.Series) -> pd.Series:
    """Percentile rank in [0, 100] (higher value -> higher percentile)."""
    x = pd.to_numeric(series, errors="coerce")
    return (x.rank(pct=True, method="average") * 100).round(2)


def _fmt_pct_raw(raw: pd.Series, pct: pd.Series) -> pd.Series:
    out = []
    for r, p in zip(raw, pct):
        if pd.isna(r):
            out.append("")
        elif pd.isna(p):
            out.append(str(r))
        else:
            out.append(f"{int(round(float(p)))}% ({r})")
    return pd.Series(out, index=raw.index, dtype="object")


def _mismatch_buckets(mm_counts: pd.Series) -> pd.DataFrame:
    """Split the ';'-joined Mode-1 mismatch_counts into mm0/mm1/mm2 columns
    (enough for the itemRgb hierarchy). Missing/blank -> 0."""

    def _get(cell, i):
        if not isinstance(cell, str) or not cell:
            return 0
        parts = cell.split(";")
        return int(parts[i]) if i < len(parts) and parts[i].isdigit() else 0

    return pd.DataFrame(
        {
            "mm0": mm_counts.map(lambda c: _get(c, 0)),
            "mm1": mm_counts.map(lambda c: _get(c, 1)),
            "mm2": mm_counts.map(lambda c: _get(c, 2)),
        }
    )


def _item_rgb(dropped, mm0, mm1, mm2, spec_pct, enpam_pct) -> np.ndarray:
    """Color hierarchy (mirrors parasol_scripts/combine_on_off_tsvs.choose_itemRgb):
    grey for non-unique (0-mm), dark grey for close (1/2-mm) off-targets, dark for
    low specificity, else a red/yellow/green ramp on the enPAM+GB percentile."""
    n = len(dropped)
    out = np.full(n, "230,106,113", dtype=object)  # default red
    enp = np.where(np.isnan(enpam_pct), 100.0, enpam_pct)
    out = np.where(enp <= 50, "230,106,113", np.where(enp <= 75, "251,243,34", "105,181,20"))
    out = np.where((~np.isnan(spec_pct)) & (spec_pct <= 55), "81,81,80", out)
    out = np.where((mm1 > 0) | (mm2 > 0), "120,120,120", out)
    out = np.where(dropped | (mm0 > 0), "150,150,150", out)
    return out


def _have_tool(name: str) -> Optional[str]:
    return shutil.which(name)


# ----------------------------------------------------------------------------
# crisprDetails.tab (+ byte offsets)
# ----------------------------------------------------------------------------
def _build_details(
    mode2_path: str,
    sidecar_path: str,
    mode1: pd.DataFrame,
    details_path: str,
    list_cap: int,
    blank_threshold: int,
) -> Dict[str, int]:
    """Write crisprDetails.tab and return ``id -> uncompressed byte offset``.

    **Streams** an arbitrarily large Mode-2 TSV: external-sort by guide_id, then
    group on the fly, so peak memory is one guide's off-targets — not the whole
    (tens-of-GB) file. One line per non-dropped guide with off_target_count > 0:
    ``mismatchCounts<TAB>offtargetList`` (browser ``pos``/strand/round(cfd*1000),
    score desc). ``list_cap <= 0`` keeps the full list; ``blank_threshold <= 0``
    never blanks. Guides with no line get _offset 0 ("No off-targets found").
    """
    import shlex

    side = pd.read_csv(sidecar_path, sep="\t", usecols=["guide_id", "id"])
    gid_to_id = dict(zip(side["guide_id"].to_numpy().tolist(), side["id"].to_numpy().tolist()))

    m1 = mode1.set_index("id")
    counts_by_id = m1["mismatch_counts"].fillna("").astype(str).to_dict()
    oc_by_id = pd.to_numeric(m1["off_target_count"], errors="coerce").to_dict()

    cap = list_cap if list_cap and list_cap > 0 else None
    blank = blank_threshold if blank_threshold and blank_threshold > 0 else None
    offsets: Dict[str, int] = {}
    seen: set = set()
    have_ot = os.path.exists(mode2_path) and os.path.getsize(mode2_path) > 0

    with open(details_path, "wb") as out:
        out.write(b"_mismatchCounts\t_crisprOfftargets\n")

        def _emit(gid_str: str, lst: str) -> None:
            counts = counts_by_id.get(gid_str, "").replace(";", ",") or "0"
            offsets[gid_str] = out.tell()
            out.write(f"{counts}\t{lst}\n".encode("ascii"))

        if have_ot:
            # External-sort the Mode-2 TSV by guide_id (numeric first field;
            # default whitespace split works since chrom names have no spaces).
            # Temp files land on the output fs; then stream + group one guide at
            # a time so memory stays bounded regardless of file size.
            sorted_path = mode2_path + ".byguide"
            tmp = os.path.dirname(os.path.abspath(details_path)) or "."
            cmd = (
                f"set -o pipefail; tail -n +2 {shlex.quote(mode2_path)} | "
                f"LC_ALL=C sort -k1,1n -S 16G -T {shlex.quote(tmp)} > {shlex.quote(sorted_path)}"
            )
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")

            def _flush(gid: int, entries: list) -> None:
                gid_str = gid_to_id.get(gid)
                if gid_str is None:
                    return
                seen.add(gid_str)
                oc = oc_by_id.get(gid_str)
                if blank is not None and oc is not None and oc > blank:
                    _emit(gid_str, "")
                    return
                entries.sort(key=lambda e: e[1], reverse=True)
                if cap is not None:
                    entries = entries[:cap]
                _emit(gid_str, "|".join(e[0] for e in entries))

            cur = None
            entries: list = []
            with open(sorted_path) as f:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    gid = int(p[0])
                    if gid != cur:
                        if cur is not None:
                            _flush(cur, entries)
                        cur, entries = gid, []
                    start, strand = int(p[2]), p[3]
                    pos = start if strand == "+" else start + 3
                    si = int(round(float(p[5]) * 1000))
                    entries.append((f"{p[1]};{pos}{strand};{si}", si))
            if cur is not None:
                _flush(cur, entries)
            os.remove(sorted_path)

        # Guides with off-targets but none in Mode-2 (all below the CFD floor):
        # keep the counts, empty list. Dropped guides (off_target_count NaN) get
        # no line -> _offset 0.
        for gid_str, oc in oc_by_id.items():
            if gid_str in seen or oc is None or pd.isna(oc) or oc <= 0:
                continue
            _emit(gid_str, "")
    return offsets


# ----------------------------------------------------------------------------
# main entry
# ----------------------------------------------------------------------------
def build_track(
    guide_df: pd.DataFrame,
    enum_prefix: str,
    outdir: str,
    chrom_sizes: Optional[str] = None,
    list_cap: int = 100,
    blank_threshold: int = 2000,
    run_tools: bool = True,
) -> Dict[str, str]:
    """Assemble the UCSC track from ``guide_df`` + crispr-ots output at
    ``enum_prefix`` (Mode-1 = ``enum_prefix``, Mode-2 = ``enum_prefix.ot.tsv``,
    sidecar = ``enum_prefix.guides.tsv``).

    ``guide_df`` must have: ``id, chrom, start, stop, strand, guideSeq, pam`` and
    any of ``DISPLAY_SCORE_COLS``. ``start``/``stop`` bound the 23-nt spacer
    (0-based, + strand). Returns a dict of the artifacts written.
    """
    os.makedirs(outdir, exist_ok=True)
    mode1 = pd.read_csv(enum_prefix)
    sidecar_path = enum_prefix + ".guides.tsv"
    mode2_path = enum_prefix + ".ot.tsv"

    score_cols = [c for c in DISPLAY_SCORE_COLS if c in guide_df.columns]

    # --- crisprDetails.tab + offsets ---
    details_path = os.path.join(outdir, DETAILS_NAME)
    offsets = _build_details(mode2_path, sidecar_path, mode1, details_path, list_cap, blank_threshold)

    # --- join guides with Mode-1 aggregates on the composite id ---
    m1 = mode1.set_index("id")
    g = guide_df.copy()
    for col in ["specificity_tttv", "specificity_tttn", "off_target_count", "mismatch_counts", "dropped"]:
        g[col] = g["id"].map(m1[col]) if col in m1.columns else np.nan
    dropped = g["dropped"].fillna(0).to_numpy() == 1

    # --- coordinates: spacer BED -> 27-nt site (5'-PAM, strand-aware) ---
    start = pd.to_numeric(g["start"], errors="coerce").fillna(0).astype(int).to_numpy()
    stop = pd.to_numeric(g["stop"], errors="coerce").fillna(0).astype(int).to_numpy()
    is_plus = g["strand"].astype(str).to_numpy() == "+"
    chrom_start = np.maximum(np.where(is_plus, start - PAM_LEN, start), 0)
    chrom_end = np.where(is_plus, stop, stop + PAM_LEN)
    thick_start = np.where(is_plus, chrom_start + PAM_LEN, chrom_start)
    thick_end = np.where(is_plus, chrom_end, chrom_end - PAM_LEN)

    # --- scores / colors / display strings ---
    # Operative specificity = TTTN: guides are TTTV-only but off-targets are
    # scored against TTTT too (enCas12a cleaves TTTT), so the bigBed score/color
    # come from specificity_tttn (Sigma cfd over ALL off-targets incl TTTT). Both
    # TTTN and TTTV are still emitted as display columns below.
    spec_tttv = pd.to_numeric(g.get("specificity_tttv"), errors="coerce")
    spec_tttn = pd.to_numeric(g.get("specificity_tttn"), errors="coerce")
    spec_pct = _pct(spec_tttn).to_numpy()
    score = np.clip(np.rint(spec_tttn.fillna(0).to_numpy() * 100), 0, 1000).astype(int)

    buckets = _mismatch_buckets(g["mismatch_counts"])
    enpam_pct = _pct(g["enpam_gb_score"]).to_numpy() if "enpam_gb_score" in g.columns else np.full(len(g), np.nan)
    item_rgb = _item_rgb(
        dropped, buckets["mm0"].to_numpy(), buckets["mm1"].to_numpy(), buckets["mm2"].to_numpy(), spec_pct, enpam_pct
    )

    # unique flags: non-dropped guides cleared the --threshold 0 screen (no 0-mm
    # off-target) -> unique; dropped (non-unique) -> False. (TTTV vs TTTN nuance
    # for dropped guides is not distinguished — they are skipped at off-target.)
    unique = np.where(dropped, "False", "True")

    out = {
        "#chr": g["chrom"].to_numpy(),
        "start": chrom_start,
        "stop": chrom_end,
        "name": " ",
        "score": score,
        "strand": g["strand"].to_numpy(),
        "thickStart": thick_start,
        "thickEnd": thick_end,
        "itemRgb": item_rgb,
        "guideSeq": g["guideSeq"].to_numpy(),
        "pam": g["pam"].to_numpy(),
    }
    for c in score_cols:
        out[c] = _fmt_pct_raw(g[c], _pct(g[c])).to_numpy()
    out["TTTV_enCas12a_specificity"] = _fmt_pct_raw(spec_tttv, _pct(spec_tttv)).to_numpy()
    out["TTTN_enCas12a_specificity"] = _fmt_pct_raw(spec_tttn, _pct(spec_tttn)).to_numpy()
    out["unique_TTTV"] = unique
    out["unique_TTTN"] = unique

    dcpf_pct = _pct(g["deepcpf1_score"]).to_numpy() if "deepcpf1_score" in g.columns else np.full(len(g), np.nan)
    mouse = [
        f"EnCas12A Spec: {'' if np.isnan(spec_pct[i]) else int(round(spec_pct[i]))}, "
        f"DeepCpf1: {'' if np.isnan(dcpf_pct[i]) else int(round(dcpf_pct[i]))}, "
        f"EnCas12a: {'' if np.isnan(enpam_pct[i]) else int(round(enpam_pct[i]))}"
        for i in range(len(g))
    ]
    out["_mouseOver"] = mouse
    out["_offset"] = [offsets.get(gid, 0) for gid in g["id"].to_numpy()]

    bed = pd.DataFrame(out)

    # --- write BED (sorted), autoSql ---
    bed_path = os.path.join(outdir, "guides.bed")
    bed_sorted = bed.sort_values(["#chr", "start"])
    bed_sorted.to_csv(bed_path, sep="\t", index=False, header=False)
    as_path = os.path.join(outdir, AS_NAME)
    with open(as_path, "w") as fh:
        fh.write(build_autosql(score_cols))

    artifacts = {"bed": bed_path, "details_tab": details_path, "autosql": as_path}

    # --- run UCSC tools if available ---
    if run_tools:
        bgzip = _have_tool("bgzip")
        if bgzip:
            subprocess.run([bgzip, "-f", "-i", details_path], check=True)
            artifacts["details_gz"] = details_path + ".gz"
            artifacts["details_gzi"] = details_path + ".gz.gzi"
        bedtobigbed = _have_tool("bedToBigBed")
        if bedtobigbed and chrom_sizes:
            bb_path = os.path.join(outdir, BB_NAME)
            subprocess.run(
                [bedtobigbed, "-tab", "-type=bed9+", f"-as={as_path}", bed_path, chrom_sizes, bb_path],
                check=True,
            )
            artifacts["bigbed"] = bb_path
    return artifacts


def _cli() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="Build the UCSC Cas12a track from crispr-ots output.")
    ap.add_argument("--guides", required=True, help="TSV with id,chrom,start,stop,strand,guideSeq,pam + score cols")
    ap.add_argument("--enum-prefix", required=True, help="crispr-ots Mode-1 output path (Mode-2/sidecar inferred)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--chrom-sizes", default=None)
    ap.add_argument("--list-cap", type=int, default=100)
    ap.add_argument("--blank-threshold", type=int, default=2000)
    ap.add_argument("--no-tools", action="store_true", help="Skip bgzip/bedToBigBed (just write .bed/.tab/.as)")
    a = ap.parse_args()
    gdf = pd.read_csv(a.guides, sep="\t", dtype=str)
    art = build_track(gdf, a.enum_prefix, a.outdir, a.chrom_sizes, a.list_cap, a.blank_threshold, not a.no_tools)
    for k, v in art.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    _cli()
