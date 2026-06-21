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
* Off-target ``pos`` = the 0-based leftmost base of the full **27-nt** site
  (4-nt PAM + 23-nt protospacer) on the forward strand = crispr-ots's ``start``,
  for **both strands**. The viewer fetches ``[pos, pos+27)``, reverse-complements
  on the ``-`` strand, and reads a 4-nt PAM + the full 23-nt match.

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

# On-target display score columns, in track column order (enseq, then deepcpf1,
# then enpam_gb). Any subset actually present in the guide table is emitted.
DISPLAY_SCORE_COLS = ["ascas12a_deepcpf1_score", "enseq_deepcpf1_score", "deepcpf1_score", "enpam_gb_score"]

# Track filenames (must not be renamed — the browser config references them).
BB_NAME = "minimumCas12A.bb"
DETAILS_NAME = "crisprDetails.tab"
AS_NAME = "cas12aTargets.as"


# ----------------------------------------------------------------------------
# autoSql
# ----------------------------------------------------------------------------
def build_autosql(
    score_cols: Sequence[str],
    spec_matrices: Sequence[tuple[str, str]] = (("enCas12a", "EnCas12a"),),
) -> str:
    """autoSql for the bigBed: bed9 + guideSeq/pam + unique flags + per-matrix
    TTTV/TTTN specificity + on-target efficiency columns + the required
    ``_mouseOver``/``_offset`` trailer. ``spec_matrices`` is ``(field_label,
    display_name)`` per off-target matrix (one TTTV + one TTTN column each). The
    column order here MUST match the BED written by :func:`build_track`."""
    score_descs = {
        "enseq_deepcpf1_score": "EnCas12a-DeepCpf1",
        "deepcpf1_score": "DeepCpf1",
        "enpam_gb_score": "EnPAM-GB",
        "ascas12a_deepcpf1_score": "AsCas12a-DeepCpf1",
    }
    head = [
        "table Cas12aTargets",
        '"CRISPR Cas12a guides, genome wide (bed9 + extra; off-targets in external bgzip tab file)"',
        "    (",
        '    string  chrom;        "Reference sequence chromosome or scaffold"',
        '    uint    chromStart;   "0-based start of the 27nt PAM+spacer site"',
        '    uint    chromEnd;     "chromStart + 27"',
        '    string  name;         "Guide ID (may be empty)"',
        '    uint    score;        "Score 0-1000 (labeled on the details page by trackDb scoreLabel)"',
        '    char[1] strand;       "+ or -"',
        '    uint    thickStart;   "Start of protospacer (thick);"',
        '    uint    thickEnd;     "End of protospacer (thick)"',
        '    uint    itemRgb;      "Display color R,G,B"',
        '    string  guideSeq;     "Guide sequence"',
        '    string  pam;          "Protospacer Adjacent Motif (PAM)"',
    ]
    # Extra string fields: unique flags, per-matrix TTTV/TTTN specificity, then
    # the on-target efficiencies. Descriptions align to the widest field name.
    extra = [("unique_TTTV", "TTTV PAM unique"), ("unique_TTTN", "TTTN PAM unique")]
    for field, display in spec_matrices:
        extra.append((f"TTTV_{field}_specificity", f"{display} TTTV CFD"))
        extra.append((f"TTTN_{field}_specificity", f"{display} TTTN CFD"))
    for c in score_cols:
        extra.append((c, score_descs.get(c, "efficiency")))
    width = max(len(f) for f, _ in extra) + 2  # "field;" + >=1 space before the description
    body = [f'    string  {(f + ";").ljust(width)}"{d}"' for f, d in extra]
    tail = [
        '    string  _mouseOver;   "Hover label (mouseOver)"',
        '    bigint  _offset;      "Uncompressed byte offset into crisprDetails.tab, or 0 if non-unique"',
        "    )",
        "",
    ]
    return "\n".join(head + body + tail)


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def _pct(series: pd.Series) -> pd.Series:
    """Percentile rank in [0, 100] (higher value -> higher percentile)."""
    x = pd.to_numeric(series, errors="coerce")
    return (x.rank(pct=True, method="average") * 100).round(2)


def _fmt_pct_raw(raw: pd.Series, pct: pd.Series) -> pd.Series:
    def _r4(v):
        """Raw value rounded to <=4 decimals; pass non-numeric through unchanged."""
        try:
            return f"{round(float(v), 4)}"
        except (TypeError, ValueError):
            return str(v)

    out = []
    for r, p in zip(raw, pct):
        if pd.isna(r):
            out.append("")
        elif pd.isna(p):
            out.append(_r4(r))
        else:
            out.append(f"{int(round(float(p)))}% ({_r4(r)})")
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


def _item_rgb(dropped, mm0, mm1, mm2, spec_pct, eff_pct) -> np.ndarray:
    """Color hierarchy (mirrors parasol_scripts/combine_on_off_tsvs.choose_itemRgb):
    grey for non-unique (0-mm), dark grey for close (1/2-mm) off-targets, dark for
    low specificity, else a red/yellow/green ramp on the AsCas12a-DeepCpf1 percentile."""
    n = len(dropped)
    out = np.full(n, "255,100,100", dtype=object)  # default red
    enp = np.where(np.isnan(eff_pct), 100.0, eff_pct)
    out = np.where(enp <= 50, "255,100,100", np.where(enp <= 75, "255,255,0", "0,200,0"))
    out = np.where((~np.isnan(spec_pct)) & (spec_pct <= 55), "80,80,80", out)
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
                    # pos = 0-based leftmost base of the full 27-nt site (4-nt PAM
                    # + 23-nt protospacer) on the forward strand = crispr-ots's
                    # `start`, for BOTH strands. The viewer fetches [pos, pos+27),
                    # reverse-complements if '-', and reads 4-nt PAM + 23-nt match.
                    start, strand = int(p[2]), p[3]
                    pos = start
                    si = int(round(float(p[5]) * 1000))
                    # Single-pass dual scoring: the Mode-2 TSV has a 2nd CFD column
                    # (cfd2). Per the spec, each off-target is
                    # chrom;posStrand;scoreInt1[;scoreInt2] (EnCas12a, then WT/2xNLS).
                    if len(p) > 6 and p[6] != "":
                        si2 = int(round(float(p[6]) * 1000))
                        entry = f"{p[1]};{pos}{strand};{si};{si2}"
                    else:
                        entry = f"{p[1]};{pos}{strand};{si}"
                    entries.append((entry, si))
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

    # Single-pass dual scoring: the secondary matrix's specificity rides in the
    # SAME Mode-1 (columns specificity_tttn_2xnls / specificity_tttv_2xnls). Emit
    # the 2xNLS-Cas12a display columns when present.
    extra_spec_data = []  # (field_label, display_name, spec_tttv, spec_tttn)
    if "specificity_tttn_2xnls" in m1.columns or "specificity_tttv_2xnls" in m1.columns:
        ids = g["id"]
        etttv = (
            pd.to_numeric(ids.map(m1["specificity_tttv_2xnls"]), errors="coerce")
            if "specificity_tttv_2xnls" in m1.columns
            else pd.Series(np.nan, index=g.index)
        )
        etttn = (
            pd.to_numeric(ids.map(m1["specificity_tttn_2xnls"]), errors="coerce")
            if "specificity_tttn_2xnls" in m1.columns
            else pd.Series(np.nan, index=g.index)
        )
        extra_spec_data.append(("AsCas12a", "AsCas12a", etttv, etttn))

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
    # Column 5 (score) + the UCSC score filter use AsCas12a (WT) TTTN specificity
    # ×100 when the secondary matrix is present, else EnCas12a TTTN. extra_spec_data
    # rows are (field, display, spec_tttv, spec_tttn) -> [0][3] is the AsCas12a TTTN.
    score_spec_tttn = extra_spec_data[0][3] if extra_spec_data else spec_tttn
    score = np.clip(np.rint(score_spec_tttn.fillna(0).to_numpy() * 100), 0, 1000).astype(int)
    # The low-specificity greying band keys off the SAME (AsCas12a) spec as the score
    # column; the hover still reports EnCas12a Spec separately via spec_pct.
    color_spec_pct = _pct(score_spec_tttn).to_numpy()

    buckets = _mismatch_buckets(g["mismatch_counts"])
    enpam_pct = _pct(g["enpam_gb_score"]).to_numpy() if "enpam_gb_score" in g.columns else np.full(len(g), np.nan)
    enseq_pct = (
        _pct(g["enseq_deepcpf1_score"]).to_numpy() if "enseq_deepcpf1_score" in g.columns else np.full(len(g), np.nan)
    )
    has_ascas12a = "ascas12a_deepcpf1_score" in g.columns
    ascas12a_pct = _pct(g["ascas12a_deepcpf1_score"]).to_numpy() if has_ascas12a else np.full(len(g), np.nan)
    # itemRgb efficiency ramp is driven by AsCas12a-DeepCpf1 (falls back to
    # EnCas12a-DeepCpf1 when the AsCas12a score isn't present). Also shown on hover.
    color_pct = ascas12a_pct if has_ascas12a else enseq_pct
    item_rgb = _item_rgb(
        dropped,
        buckets["mm0"].to_numpy(),
        buckets["mm1"].to_numpy(),
        buckets["mm2"].to_numpy(),
        color_spec_pct,
        color_pct,
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
    # Column order must match build_autosql: unique flags, then TTTV/TTTN
    # specificity, then the on-target efficiency scores (enseq, deepcpf1, enpam_gb).
    out["unique_TTTV"] = unique
    out["unique_TTTN"] = unique
    out["TTTV_enCas12a_specificity"] = _fmt_pct_raw(spec_tttv, _pct(spec_tttv)).to_numpy()
    out["TTTN_enCas12a_specificity"] = _fmt_pct_raw(spec_tttn, _pct(spec_tttn)).to_numpy()
    for field, _display, etttv, etttn in extra_spec_data:
        out[f"TTTV_{field}_specificity"] = _fmt_pct_raw(etttv, _pct(etttv)).to_numpy()
        out[f"TTTN_{field}_specificity"] = _fmt_pct_raw(etttn, _pct(etttn)).to_numpy()
    for c in score_cols:
        out[c] = _fmt_pct_raw(g[c], _pct(g[c])).to_numpy()

    dcpf_pct = _pct(g["deepcpf1_score"]).to_numpy() if "deepcpf1_score" in g.columns else np.full(len(g), np.nan)

    def _pf(v):
        return "" if np.isnan(v) else int(round(v))

    # percentiles of each extra matrix's TTTN specificity (for the hover)
    extra_spec_pct = [(display, _pct(etttn).to_numpy()) for _field, display, _etttv, etttn in extra_spec_data]

    def _hover(i):
        # operative (primary) specificity; non-unique guides have none.
        parts = ["Not unique at TTTN" if np.isnan(spec_pct[i]) else f"EnCas12a Spec: {int(round(spec_pct[i]))}"]
        for display, pcts in extra_spec_pct:
            if not np.isnan(pcts[i]):
                parts.append(f"{display} Spec: {int(round(pcts[i]))}")
        parts += [
            f"EnCas12a-DeepCpf1: {_pf(enseq_pct[i])}",
            f"DeepCpf1: {_pf(dcpf_pct[i])}",
            f"EnPAM-GB: {_pf(enpam_pct[i])}",
        ]
        if has_ascas12a:
            parts.append(f"AsCas12a-DeepCpf1: {_pf(ascas12a_pct[i])}")
        return ", ".join(parts)

    mouse = [_hover(i) for i in range(len(g))]
    out["_mouseOver"] = mouse
    out["_offset"] = [offsets.get(gid, 0) for gid in g["id"].to_numpy()]

    bed = pd.DataFrame(out)

    # --- write BED (sorted), autoSql ---
    bed_path = os.path.join(outdir, "guides.bed")
    bed_sorted = bed.sort_values(["#chr", "start"])
    bed_sorted.to_csv(bed_path, sep="\t", index=False, header=False)
    spec_matrices = [("enCas12a", "EnCas12a")] + [(field, display) for field, display, _tv, _tn in extra_spec_data]
    as_path = os.path.join(outdir, AS_NAME)
    with open(as_path, "w") as fh:
        fh.write(build_autosql(score_cols, spec_matrices))

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


def merge_ucscgb_tracks(
    track_dirs: Sequence[str],
    out: str,
    chrom_sizes: str,
    as_file: str,
    bedtobigbed: str = "bedToBigBed",
    sort_mem: str = "8G",
    sort_threads: str = "8",
    keep_bed: bool = False,
) -> Dict[str, str]:
    """Merge per-chunk Cas12a UCSC tracks (each produced by a single-pass
    :func:`build_track`) into one track.

    Each ``track_dir`` must contain ``guides.bed`` (bed9+, last column =
    ``_offset`` = uncompressed byte offset into that piece's
    ``crisprDetails.tab``; ``0`` = non-unique sentinel) and
    ``crisprDetails.tab.gz``. The off-target details carry no guide key, so the
    BED ``_offset`` is the only link: we stream each piece's details, map every
    data line's old byte offset to its offset in the merged file, rewrite that
    piece's BED ``_offset``, then sort + ``bedToBigBed`` and ``bgzip`` + index.
    One shared header sits at byte 0 (the ``0`` sentinel). Merge order is
    arbitrary -- offsets are reassigned here and the bigBed is sorted at the end.
    Returns the artifact paths.
    """
    header = b"_mismatchCounts\t_crisprOfftargets\n"

    def _zcat_lines(path):
        p = subprocess.Popen(["bgzip", "-dc", "-@", "4", path], stdout=subprocess.PIPE, bufsize=1024 * 1024)
        try:
            for line in p.stdout:
                yield line
        finally:
            p.stdout.close()
            if p.wait() != 0:
                raise RuntimeError(f"bgzip -dc failed for {path}")

    os.makedirs(out, exist_ok=True)
    merged_tab = os.path.join(out, DETAILS_NAME)
    merged_bed = os.path.join(out, "guides.merged.bed")

    n_guides = n_unique = 0
    running = 0  # byte offset in the merged uncompressed details
    with open(merged_tab, "wb") as outf, open(merged_bed, "wb") as bedout:
        outf.write(header)
        running = len(header)
        for tdir in track_dirs:
            det = os.path.join(tdir, DETAILS_NAME + ".gz")
            bed = os.path.join(tdir, "guides.bed")
            if not (os.path.exists(det) and os.path.exists(bed)):
                raise FileNotFoundError(f"missing track files in {tdir}")
            # 1. stream this piece's details -> old(byte in piece) -> new(byte in merged)
            offmap = {}
            piece_off = 0
            for i, line in enumerate(_zcat_lines(det)):
                if i == 0:  # per-piece header at byte 0 -> not copied (merged has its own at 0)
                    if line != header:
                        raise ValueError(f"unexpected details header in {det}: {line!r}")
                    piece_off += len(line)
                    continue
                offmap[piece_off] = running
                outf.write(line)
                running += len(line)
                piece_off += len(line)
            # 2. rewrite this piece's BED with remapped _offset (last column)
            with open(bed, "rb") as bf:
                for row in bf:
                    row = row.rstrip(b"\n")
                    if not row:
                        continue
                    fields = row.split(b"\t")
                    old = int(fields[-1])
                    if old == 0:
                        new = 0
                    else:
                        try:
                            new = offmap[old]
                        except KeyError:
                            raise KeyError(f"BED _offset {old} not found in {det} (offset map mismatch)")
                        n_unique += 1
                    fields[-1] = str(new).encode()
                    bedout.write(b"\t".join(fields) + b"\n")
                    n_guides += 1
            offmap.clear()
            print(f"  merged {os.path.basename(tdir)}: cum_guides={n_guides} merged_bytes={running}", flush=True)

    # 3. sort the merged BED and build the bigBed
    sorted_bed = os.path.join(out, "guides.sorted.bed")
    subprocess.run(
        f"LC_ALL=C sort -S{sort_mem} --parallel={sort_threads} -T {out} -k1,1 -k2,2n {merged_bed} > {sorted_bed}",
        shell=True,
        check=True,
    )
    bb = os.path.join(out, BB_NAME)
    subprocess.run([bedtobigbed, "-tab", "-type=bed9+", f"-as={as_file}", sorted_bed, chrom_sizes, bb], check=True)
    # 4. bgzip + index the merged details (offsets are uncompressed; .gzi handles the seek)
    subprocess.run(["bgzip", "-f", "-@", "8", merged_tab], check=True)
    subprocess.run(["bgzip", "-r", merged_tab + ".gz"], check=True)
    shutil.copy(as_file, os.path.join(out, AS_NAME))
    os.remove(merged_bed)
    if not keep_bed:
        os.remove(sorted_bed)
    print(
        f"\nMERGE DONE -> {out}\n  guides={n_guides}  unique(with details)={n_unique}  details_bytes={running}",
        flush=True,
    )
    return {
        "bigbed": bb,
        "details_gz": merged_tab + ".gz",
        "details_gzi": merged_tab + ".gz.gzi",
        "autosql": os.path.join(out, AS_NAME),
    }


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
