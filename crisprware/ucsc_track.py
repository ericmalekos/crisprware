#!/usr/bin/env python3
"""Build a UCSC Genome Browser Cas12a track from crispr-ots scored output.

Consumes the per-guide on-target scores plus the crispr-ots streaming
off-target output (Mode-1 aggregated CSV + Mode-2 per-off-target TSV + the
guide sidecar) and emits the three files pinned by
``cas12a_track_output_spec.md``:

* ``cas12a.bb``             bigBed (bed9+) of on-target guides
* ``crisprDetails.tab.gz``  bgzip per-guide off-target details (+ ``.gzi``)
* ``cas12aTargets.as``      autoSql schema (build input)

Coordinate conventions (validated against hg38, both strands):

* On-target feature spans PAM(4) + spacer(23) = 27 nt. The guide BED's
  ``start``/``stop`` bound the 23-nt spacer on the + strand; the 4-nt PAM is
  added on the 5' side per strand (``+``: left, ``-``: right).
* Off-target ``pos`` = the 0-based leftmost base of the full **27-nt** site
  (4-nt PAM + 23-nt protospacer) on the forward strand = crispr-ots's ``start``,
  for **both strands**. The viewer fetches ``[pos, pos+27)``, reverse-complements
  on the ``-`` strand, and reads a 4-nt PAM + the full 23-nt match.

``_mismatchCounts`` and ``_crisprOfftargets`` describe the SAME population: both
are built from the Mode-2 stream (floored at the run's ``--cfd-threshold``), so
``sum(_mismatchCounts) == len(_crisprOfftargets)`` for every listed guide. This
is the Cas9/CRISPOR convention: an off-target below the floor is neither listed
nor counted. The exception, also per Cas9, is a guide with too many off-targets
to show (over the per-guide cap, or over ``blank_threshold``): the list is
emitted empty and the histogram still reports the true above-floor totals, which
is what makes the browser print its "too many off-targets" message.

Per-guide specificity is deliberately NOT reconstructible from these two fields:
it sums every off-target including sub-floor ones (as MIT specificity does for
Cas9).
"""

from __future__ import annotations

import csv
import mmap
import os
import shutil
import subprocess
from importlib import resources
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

SITE_LEN = 27
SPACER_LEN = 23
PAM_LEN = 4


# ----------------------------------------------------------------------------
# Genome PAM reader (off-target TTTV/TTTN classification, shared with the
# re-derive script). AsCas12a does not cleave TTTT-prefixed off-targets, so the
# AsCas12a per-off-target CFD is zeroed for them (Layer 2) and they are excluded
# from the AsCas12a close-off-target grey (Layer 3). The scanner is PAM-blind for
# the per-off-target CFD, so the PAM is recovered here from the reference FASTA.
# ----------------------------------------------------------------------------
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


class Genome:
    """Minimal mmap FASTA reader (needs a ``.fai``). Output is UPPERCASED so
    soft-masked (lowercase) repeats -- exactly where TTTT PAMs concentrate -- are
    classified correctly. Random access is O(1); safe for genome-scale streaming."""

    def __init__(self, fna: str):
        self._fh = open(fna, "rb")
        self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)
        self._idx: dict = {}
        with open(fna + ".fai") as f:
            for line in f:
                name, length, offset, lb, lw = line.rstrip("\n").split("\t")[:5]
                self._idx[name] = (int(offset), int(lb), int(lw), int(length))

    def _byte(self, meta, p: int) -> int:
        offset, lb, lw, _ = meta
        return offset + (p // lb) * lw + (p % lb)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """0-based half-open [start, end) on the + strand, uppercase, newline-stripped."""
        meta = self._idx.get(chrom)
        if meta is None:
            return ""
        start, end = max(0, start), min(meta[3], end)
        if start >= end:
            return ""
        b, e = self._byte(meta, start), self._byte(meta, end - 1) + 1
        return self._mm[b:e].replace(b"\n", b"").decode("ascii", "replace").upper()

    def pam(self, chrom: str, pos: int, strand: str) -> str:
        """4-nt 5' Cas12a PAM of an off-target whose 27-nt site starts at ``pos``
        (``+``: [pos, pos+4); ``-``: revcomp of [pos+23, pos+27))."""
        if strand == "+":
            return self.fetch(chrom, pos, pos + 4)
        return _revcomp(self.fetch(chrom, pos + 23, pos + 27))

    def spacer(self, chrom: str, pos: int, strand: str) -> str:
        """23-nt protospacer of an off-target whose 27-nt site starts at ``pos``
        (``+``: [pos+4, pos+27); ``-``: revcomp of [pos, pos+23)). Same geometry as
        the on-target site; reproduces the guide's stored guideSeq. The pipeline
        drives Layer 3 off the scanner's `mm` field, so this is not used at
        generation, but the re-derive imports it to re-Hamming from the genome."""
        if strand == "+":
            return self.fetch(chrom, pos + 4, pos + 27)
        return _revcomp(self.fetch(chrom, pos, pos + 23))


def is_tttt(pam: str) -> bool:
    """TTTT PAM -- AsCas12a does not cleave it (enCas12a does). Anything not a clean
    4-nt TTTT (Ns, or a truncated read near a contig edge) is treated as keep."""
    return pam == "TTTT"


# ----------------------------------------------------------------------------
# enCas12a off-target PAM-activity weights (for the enCas12a track's four
# guidescan specificities). The enumerate scan uses the fa>0.3 (30-PAM) index,
# so every off-target sits at a fa>0.3 PAM; this table gives, per PAM, its
# frac-active weight (relative cleavage vs the best PAM) and its fa>0.6/fa>0.3
# activity-cutoff membership.
# ----------------------------------------------------------------------------
ENCAS12A_WEIGHTS_FILE = "encas12a_offtarget_weights.tsv"


def default_encas12a_weights_path() -> str:
    """Path to the shipped enCas12a off-target PAM-activity weights table."""
    try:
        return str(resources.files("crisprware").joinpath("data", ENCAS12A_WEIGHTS_FILE))
    except (TypeError, ModuleNotFoundError, AttributeError):
        return os.path.join(os.path.dirname(os.path.dirname(__file__)), "crisprware", "data", ENCAS12A_WEIGHTS_FILE)


def load_encas12a_weights(path: Optional[str] = None) -> Dict[str, tuple]:
    """enCas12a off-target PAM weights: ``pam -> (weight, in_fa0.6, in_fa0.3)``.

    ``weight`` is the frac-active weight; ``in_fa0.6``/``in_fa0.3`` are the high/moderate
    activity-cutoff membership flags. Consumed by :func:`_build_details` to aggregate the
    four guidescan specificities (raw + PAM-weighted, at the fa>0.6 and fa>0.3 cutoffs)."""
    weights: Dict[str, tuple] = {}
    with open(path or default_encas12a_weights_path()) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            weights[r["pam"]] = (float(r["weight_fracactive"]), r["in_fa0.6"] == "1", r["in_fa0.3"] == "1")
    return weights


# On-target display score columns, in track column order (enseq, then deepcpf1,
# then enpam_gb). Any subset actually present in the guide table is emitted.
DISPLAY_SCORE_COLS = ["ascas12a_deepcpf1_score", "enseq_deepcpf1_score", "deepcpf1_score", "enpam_gb_score"]

# Track filenames; the deployed browser config must reference these same names.
BB_NAME = "cas12a.bb"
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
        extra.append((f"TTTV_{field}_specificity", f"{display} TTTV specificity"))
        extra.append((f"TTTN_{field}_specificity", f"{display} TTTN specificity"))
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


def _fmt_raw(raw: pd.Series) -> pd.Series:
    """Raw value to <=4 decimals as a string; blank for NaN. Used for the specificity
    column: a percentile rank is degenerate there because the vast majority of guides
    tie at a perfectly-specific 1.0, so the raw CFD specificity (0-1) is what's shown."""
    out = []
    for r in raw:
        out.append("" if pd.isna(r) else f"{round(float(r), 4)}")
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


# Absolute AsCas12a-TTTV specificity thresholds for the two-tier low-specificity
# greys. TTTV (not TTTN) because AsCas12a does not cleave TTTT-prefixed sites, and
# ABSOLUTE (not a percentile rank) so a genome of uniformly specific guides (e.g.
# yeast, where ~93% would otherwise be greyed despite median specificity 1.0) is
# colored by efficiency. Mirrors CRISPOR's "specificity > 50 is good"; the browser
# score filter is the same raw x100.
SPEC_GREY_BAD = 0.20  # raw AsCas12a-TTTV < 0.20      -> darkest grey (bad)
SPEC_GREY_MARGINAL = 0.50  # 0.20 <= raw < 0.50       -> lighter grey (marginal)

# Efficiency-percentile cuts for the red/yellow/green ramp, taken from the Cas9
# crisprAllTargets track description so the two tracks read the same way:
# "low predicted cleavage: ... percentile <= 30", "medium ...: > 30 and < 55",
# "high ...: > 55". (These were 50/75 before, markedly stricter than Cas9, so a
# guide at percentile 60 showed yellow here and green on the Cas9 track.)
EFF_RAMP_LOW = 30  # percentile <= 30            -> red
EFF_RAMP_HIGH = 55  # 30 < percentile <= 55      -> yellow; > 55 -> green

# Cuts for the specificity-only _specColor, read off the Cas9 crisprAllTargets sample
# rows (specificity 27/37/43 -> red, 51/52/65/69 -> olive, 72/76 -> green). On a
# 0-100 scale, so the raw 0-1 AsCas12a specificity is multiplied by 100 first.
SPEC_COLOR_BAD = 50  # <= 50      -> red
SPEC_COLOR_MID = 70  # 50 to 70   -> olive; > 70 -> green


def _cas9_style_name(col: str) -> str:
    """``deepcpf1_score`` -> ``deepcpf1``; ``ascas12a_deepcpf1_score`` ->
    ``ascas12aDeepcpf1``. Cas9 names each efficiency model with a bare camelCase
    field plus a ``...Perc`` sibling (``fusi``/``fusiPerc``,
    ``crisprScan``/``crisprScanPerc``); this derives that base name."""
    s = col[:-6] if col.endswith("_score") else col
    head, *rest = s.split("_")
    return head + "".join(p.capitalize() for p in rest)


def _spec_eff_color(spec_raw, eff_pct) -> np.ndarray:
    """Color for a unique guide with no close (1/2-mm) off-targets: two-tier
    ABSOLUTE low-specificity grey on the raw AsCas12a-TTTV specificity, else the
    red/yellow/green efficiency ramp on the AsCas12a-DeepCpf1 percentile. Missing
    specificity is treated as fully specific (not greyed)."""
    enp = np.where(np.isnan(eff_pct), 100.0, eff_pct)
    out = np.where(
        enp <= EFF_RAMP_LOW, "255,100,100", np.where(enp <= EFF_RAMP_HIGH, "255,255,0", "0,200,0")
    )
    sr = np.where(np.isnan(spec_raw), 1.0, spec_raw)
    out = np.where(sr < SPEC_GREY_MARGINAL, "90,90,90", out)  # marginal 0.20-0.50
    out = np.where(sr < SPEC_GREY_BAD, "50,50,50", out)  # bad < 0.20 (overrides marginal)
    return out


def _item_rgb(dropped, mm0, mm1, mm2, spec_raw, eff_pct) -> np.ndarray:
    """Color hierarchy (highest priority applied last): efficiency ramp / two-tier
    absolute low-specificity grey (see _spec_eff_color), then grey for close (1/2-mm)
    off-targets, then light grey for non-unique (0-mm / dropped). ``spec_raw`` is the
    raw AsCas12a-TTTV specificity; ``mm*`` are the (TTTV-restricted for AsCas12a) counts.

    Used by EVERY path including AsCas12a: this mirrors the Cas9 default scheme, where
    grey by specificity wins and the efficiency ramp shows only for specific guides."""
    out = _spec_eff_color(spec_raw, eff_pct)
    out = np.where((mm1 > 0) | (mm2 > 0), "120,120,120", out)
    out = np.where(dropped | (mm0 > 0), "150,150,150", out)
    return out


def _spec_color(dropped, mm0, mm1, mm2, spec_raw) -> np.ndarray:
    """Specificity-ONLY color, the Cas9 ``_specColor`` ALTERNATIVE coloring.

    Cas9 ships three colors per guide: ``reserved``/itemRgb (the default: grey when
    MIT specificity <= 50, otherwise a Doench-2016 efficiency-percentile ramp) plus
    two alternatives, ``_crisprScanColor`` (Moreno-Mateos) and ``_specColor`` (MIT
    specificity). The alternatives are other ways to color the SAME guide, not a
    decomposition of the default, so this column is purely additive: ``itemRgb``
    keeps the blended default and a viewer can instead color by specificity alone.

    Cas9 makes this a PURE function of specificity with no uniqueness greys: its
    sample rows pair specificity 27/37/43 with red, 51/52/65/69 with olive and 72/76
    with green, and keep red even where ``reserved`` has already greyed out. So the
    cuts are ``<= 50`` red, ``<= 70`` olive, ``> 70`` green, and a greyed guide still
    gets a colour here. ``spec_raw`` is the raw AsCas12a-TTTV specificity (0-1).
    A non-unique guide has no meaningful specificity, so it takes the worst tier;
    a missing specificity keeps the usual "treat as fully specific" convention."""
    pct = np.where(np.isnan(spec_raw), 100.0, np.asarray(spec_raw, dtype=float) * 100.0)
    pct = np.where(dropped | (mm0 > 0), 0.0, pct)  # non-unique -> worst tier
    return np.where(
        pct <= SPEC_COLOR_BAD, "255,0,0", np.where(pct <= SPEC_COLOR_MID, "128,128,0", "0,255,0")
    )


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
    genome: Optional["Genome"] = None,
    weights: Optional[Dict[str, tuple]] = None,
    oncount: Optional[Dict[str, int]] = None,
    specs_out: Optional[Dict[str, tuple]] = None,
    list_floor: float = 0.0,
):
    """Write crisprDetails.tab and return ``(offsets, tttv_mm)``.

    **enCas12a track** (``weights`` given): in the SAME Mode-2 stream, aggregate the four
    guidescan specificities per guide from each off-target's reference PAM -- raw + PAM
    weighted, at the fa>0.6 and fa>0.3 activity cutoffs -- and fill ``specs_out[id] =
    (spec_fa0.6, spec_fa0.6_w, spec_fa0.3, spec_fa0.3_w, n_ot_fa0.6, n_ot_fa0.3)`` where
    ``spec = 1 / (on_count + sum)`` (``on_count`` from ``oncount[id]``, default 1). The spec
    sums count EVERY off-target (run the scan unfloored); ``list_floor`` floors only the
    displayed popup list, keeping the specs exact. Requires ``genome`` (for the off-target
    PAM). ``weights=None`` leaves every enCas12a path inert (AsCas12a/legacy unchanged).

    **Streams** an arbitrarily large Mode-2 TSV: external-sort by guide_id, then
    group on the fly, so peak memory is one guide's off-targets, not the whole
    (tens-of-GB) file. One line per non-dropped guide with off_target_count > 0:
    ``mismatchCounts<TAB>offtargetList`` (browser ``pos``/strand/round(cfd*1000),
    score desc). ``list_cap <= 0`` keeps the full list; ``blank_threshold <= 0``
    never blanks. Guides with no line get _offset 0 ("No off-targets found").

    ``mismatchCounts`` is binned from the SAME floored off-targets that the list
    shows, not from Mode-1, so ``sum(counts) == len(list)`` per guide (Cas9
    convention; Mode-1 totals include sub-floor hits the popup would never
    display). A guide over ``list_cap`` or ``blank_threshold`` emits an empty
    list while keeping its true above-floor histogram, so the browser shows its
    "too many off-targets" message instead of a silently truncated list. Bin
    width follows Mode-1's histogram width.

    ``offsets`` maps ``id -> uncompressed byte offset``. When ``genome`` is given the
    AsCas12a per-off-target CFD (``cfd2``) is zeroed for TTTT off-targets (Layer 2;
    AsCas12a does not cleave them) and ``tttv_mm`` maps ``id -> [mm0, mm1, mm2]`` over
    the TTTV-restricted (non-TTTT) off-targets that cleared the CFD floor, for the
    AsCas12a close-off-target grey (Layer 3). Without a genome ``tttv_mm`` is empty and
    callers fall back to the Mode-1 (TTTN) mismatch counts.
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
    tttv_mm: Dict[str, list] = {}
    seen: set = set()
    have_ot = os.path.exists(mode2_path) and os.path.getsize(mode2_path) > 0

    # Histogram width = Mode-1's bin count (max_mismatches + 1); Mode-1 supplies
    # the shape only, never the totals (those are re-binned from the floored
    # Mode-2 off-targets below).
    _widths = [s.count(";") + 1 for s in counts_by_id.values() if s]
    nbins = max(_widths) if _widths else 5
    zero_hist = ",".join(["0"] * nbins)

    def _hist(mms) -> str:
        b = [0] * nbins
        for mm in mms:
            if 0 <= mm < nbins:
                b[mm] += 1
        return ",".join(str(x) for x in b)

    with open(details_path, "wb") as out:
        out.write(b"_mismatchCounts\t_crisprOfftargets\n")

        def _emit(gid_str: str, lst: str, counts: str) -> None:
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

            def _flush(gid: int, entries: list, acc: Optional[dict] = None) -> None:
                gid_str = gid_to_id.get(gid)
                if gid_str is None:
                    return
                seen.add(gid_str)
                # enCas12a: four guidescan specificities from this guide's accumulated
                # off-target CFD sums (computed regardless of popup blank/cap).
                if specs_out is not None and acc is not None:
                    oc = (oncount or {}).get(gid_str, 1)
                    specs_out[gid_str] = (
                        1.0 / (oc + acc["c06"]),
                        1.0 / (oc + acc["w06"]),
                        1.0 / (oc + acc["c03"]),
                        1.0 / (oc + acc["w03"]),
                        acc["n06"],
                        acc["n03"],
                    )
                # Layer 3: TTTV-restricted [mm0, mm1, mm2] over ALL floored
                # off-targets of this guide (before capping/blanking), TTTT excluded.
                if genome is not None:
                    b = [0, 0, 0]
                    for _e, _si, mm, tttt in entries:
                        if not tttt and 0 <= mm < 3:
                            b[mm] += 1
                    tttv_mm[gid_str] = b
                # Bin the floored off-targets this guide actually has, so the
                # histogram and the list agree (Cas9 convention).
                hist = _hist(e[2] for e in entries)
                oc = oc_by_id.get(gid_str)
                if blank is not None and oc is not None and oc > blank:
                    _emit(gid_str, "", hist)
                    return
                if cap is not None and len(entries) > cap:
                    # Too many to show. Cas9 emits an empty list here rather than a
                    # truncated one; the histogram keeps the true above-floor totals
                    # and the browser prints its "use the Crispor.org link" message.
                    _emit(gid_str, "", hist)
                    return
                entries.sort(key=lambda e: e[1], reverse=True)
                _emit(gid_str, "|".join(e[0] for e in entries), hist)

            def _new_acc() -> dict:
                return {"c06": 0.0, "w06": 0.0, "c03": 0.0, "w03": 0.0, "n06": 0, "n03": 0}

            cur = None
            entries: list = []
            acc = _new_acc() if weights is not None else None
            with open(sorted_path) as f:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    gid = int(p[0])
                    if gid != cur:
                        if cur is not None:
                            _flush(cur, entries, acc)
                        cur, entries = gid, []
                        acc = _new_acc() if weights is not None else None
                    # pos = 0-based leftmost base of the full 27-nt site (4-nt PAM
                    # + 23-nt protospacer) on the forward strand = crispr-ots's
                    # `start`, for BOTH strands. The viewer fetches [pos, pos+27),
                    # reverse-complements if '-', and reads 4-nt PAM + 23-nt match.
                    start, strand = int(p[2]), p[3]
                    pos = start
                    mm = int(p[4])
                    cfd = float(p[5])
                    si = int(round(cfd * 1000))
                    # Reference 4-nt PAM of this off-target (the scanner is PAM-blind): used
                    # for the AsCas12a TTTT check (Layer 2) and the enCas12a spec weights.
                    pam_str = genome.pam(p[1], pos, strand) if genome is not None else ""
                    tttt = is_tttt(pam_str) if genome is not None else False
                    # enCas12a: aggregate the four guidescan specificities over EVERY
                    # off-target (unfloored), keyed on the off-target PAM's activity weight
                    # + fa>0.6/fa>0.3 cutoff membership; 1/(on_count + sum) at flush.
                    if weights is not None:
                        wt = weights.get(pam_str)
                        if wt is not None:
                            w, in06, in03 = wt
                            wcfd = cfd * w
                            if in03:
                                acc["c03"] += cfd
                                acc["w03"] += wcfd
                                acc["n03"] += 1
                            if in06:
                                acc["c06"] += cfd
                                acc["w06"] += wcfd
                                acc["n06"] += 1
                    # Single-pass dual scoring: the Mode-2 TSV has a 2nd CFD column
                    # (cfd2). Per the spec, each off-target is
                    # chrom;posStrand;scoreInt1[;scoreInt2] (EnCas12a, then WT/2xNLS).
                    if len(p) > 6 and p[6] != "":
                        si2 = 0 if tttt else int(round(float(p[6]) * 1000))
                        entry = f"{p[1]};{pos}{strand};{si};{si2}"
                    else:
                        entry = f"{p[1]};{pos}{strand};{si}"
                    # list_floor floors ONLY the displayed popup (enCas12a scans unfloored so
                    # the spec sums above are exact); default 0 keeps AsCas12a/legacy behavior.
                    if cfd >= list_floor:
                        entries.append((entry, si, mm, tttt))
            if cur is not None:
                _flush(cur, entries, acc)
            os.remove(sorted_path)

        # Guides whose off-targets are ALL below the CFD floor: nothing to list and,
        # under the Cas9 convention, nothing to count either -> zero histogram, empty
        # list. Dropped guides (off_target_count NaN) get no line -> _offset 0.
        for gid_str, oc in oc_by_id.items():
            if gid_str in seen or oc is None or pd.isna(oc) or oc <= 0:
                continue
            _emit(gid_str, "", zero_hist)
    return offsets, tttv_mm


# ----------------------------------------------------------------------------
# secondary-matrix specificity columns
# ----------------------------------------------------------------------------
# Display names for known matrix-name slugs. crispr-ots derives the column
# suffix from the matrix's own name, so `2xNLS-Cas12a` -> `2xnls_cas12a`.
SECONDARY_MATRIX_LABELS = {
    "2xnls_cas12a": "AsCas12a",
    "2xnls": "AsCas12a",  # legacy: suffix was hardcoded before the rename
    "encas12a": "EnCas12a",
}


def secondary_specificity_suffixes(columns: Sequence[str]) -> list[str]:
    """Matrix-name suffixes of any secondary specificity columns, in order.

    A dual-matrix Mode-1 file carries `specificity_tttn_<matrix>` and
    `specificity_tttv_<matrix>` alongside the unsuffixed primary pair. The
    suffix used to be the hardcoded literal `2xnls` for whichever matrix was
    passed second, which mislabelled the enzyme whenever the two `--score`
    matrices were given in the other order. Matching by prefix keeps this
    working for any matrix name, and the legacy spelling still resolves.
    """
    out: list[str] = []
    for col in columns:
        for stem in ("specificity_tttv_", "specificity_tttn_"):
            if col.startswith(stem):
                suffix = col[len(stem):]
                if suffix and suffix not in out:
                    out.append(suffix)
    return out


def build_autosql_ascas12a(ontarget_cols: Sequence[str] = ()) -> str:
    """autoSql for the single-enzyme **AsCas12a** track: bed9 + ``_specColor`` +
    guideSeq/pam + a single ``unique`` flag + one ``AsCas12a_specificity`` column +
    one column per on-target model in ``ontarget_cols`` + the ``_mouseOver``/``_offset``
    trailer. AsCas12a is scored against a TTTV-only off-target index, so there is NO
    TTTN column and NO enCas12a column.

    Following the Cas9 track, the two axes get their own colors: ``itemRgb`` is the
    on-target efficiency ramp (Cas9's ``reserved``/Doench-2016 role) and ``_specColor``
    carries the specificity and uniqueness greys (Cas9's MIT-specificity role). Color
    fields sit directly after ``itemRgb``, before ``guideSeq``, as Cas9 orders them.
    Column order MUST match the BED written by :func:`build_track` with
    ``enzyme="ascas12a"``."""
    head = [
        "table AsCas12aTargets",
        '"CRISPR AsCas12a (WT/2xNLS) guides, genome wide (bed9 + extra; TTTV off-target search; off-targets in external bgzip tab file)"',
        "    (",
        '    string  chrom;        "Reference sequence chromosome or scaffold"',
        '    uint    chromStart;   "0-based start of the 27nt PAM+spacer site"',
        '    uint    chromEnd;     "chromStart + 27"',
        '    string  name;         "Guide ID (may be empty)"',
        '    uint    score;        "AsCas12a-TTTV specificity x100 (0-100)"',
        '    char[1] strand;       "+ or -"',
        '    uint    thickStart;   "Start of protospacer (thick);"',
        '    uint    thickEnd;     "End of protospacer (thick)"',
        '    uint    itemRgb;      "Display color R,G,B (specificity grey, else DeepCpf1 efficiency ramp)"',
        '    string  _crisprDeepCpf1Color; "Alternative color R,G,B by the second on-target model (seq-DeepCpf1)"',
        '    string  _specColor;   "Alternative color R,G,B by AsCas12a-TTTV specificity alone"',
        '    string  guideSeq;     "Guide sequence"',
        '    string  pam;          "Protospacer Adjacent Motif (PAM)"',
    ]
    extra = [
        ("unique", "Guide has a unique on-target (TTTV off-target search)"),
        ("AsCas12a_specificity", "AsCas12a (WT/2xNLS) CFD specificity, 0-1; TTTV off-target search"),
    ]
    ot_descs = {
        "deepcpf1_score": "DeepCpf1 (Kim 2018) AsCas12a on-target efficiency",
        "ascas12a_deepcpf1_score": "AsCas12a seq-DeepCpf1 variant (Chen 2025) on-target efficiency",
    }
    for c in ontarget_cols:
        d = ot_descs.get(c, "AsCas12a on-target efficiency")
        base = _cas9_style_name(c)
        extra.append((base, f"Efficiency: {d} (raw)"))
        extra.append((base + "Perc", f"Efficiency: {d} (percentile)"))
    width = max(len(f) for f, _ in extra) + 2
    body = [f'    string  {(f + ";").ljust(width)}"{d}"' for f, d in extra]
    tail = [
        '    string  _mouseOver;   "Hover label (mouseOver)"',
        '    bigint  _offset;      "Uncompressed byte offset into crisprDetails.tab, or 0 if non-unique"',
        "    )",
        "",
    ]
    return "\n".join(head + body + tail)


def build_autosql_encas12a(ontarget_cols: Sequence[str] = ()) -> str:
    """autoSql for the single-enzyme **enCas12a** track: bed9 + guideSeq/pam + a unique
    flag + the four guidescan specificities (raw + PAM-activity-weighted, at the fa>0.6 and
    fa>0.3 activity cutoffs) + one column per on-target model in ``ontarget_cols`` + the
    ``_mouseOver``/``_offset`` trailer. enCas12a is scored against the broad fa>0.3 (30-PAM)
    off-target index; the fa>0.3 PAM-weighted specificity drives the score/color. Column
    order MUST match the BED written by :func:`build_track` with ``enzyme="encas12a"``."""
    head = [
        "table EnCas12aTargets",
        '"CRISPR enCas12a (enAsCas12a) guides, genome wide (bed9 + extra; fa>0.3 off-target search; off-targets in external bgzip tab file)"',
        "    (",
        '    string  chrom;        "Reference sequence chromosome or scaffold"',
        '    uint    chromStart;   "0-based start of the 27nt PAM+spacer site"',
        '    uint    chromEnd;     "chromStart + 27"',
        '    string  name;         "Guide ID (may be empty)"',
        '    uint    score;        "enCas12a fa>0.3 PAM-weighted specificity x100 (0-100)"',
        '    char[1] strand;       "+ or -"',
        '    uint    thickStart;   "Start of protospacer (thick);"',
        '    uint    thickEnd;     "End of protospacer (thick)"',
        '    uint    itemRgb;      "Display color R,G,B"',
        '    string  guideSeq;     "Guide sequence"',
        '    string  pam;          "Protospacer Adjacent Motif (PAM)"',
    ]
    extra = [
        ("unique", "Guide has a unique on-target (fa>0.3 off-target search)"),
        ("enCas12a_spec_fa06", "enCas12a CFD specificity 0-1; off-targets at high-activity fa>0.6 PAMs"),
        ("enCas12a_spec_fa06_wt", "enCas12a PAM-activity-weighted specificity 0-1; fa>0.6 PAMs"),
        ("enCas12a_spec_fa03", "enCas12a CFD specificity 0-1; off-targets at moderate-activity fa>0.3 PAMs"),
        ("enCas12a_spec_fa03_wt", "enCas12a PAM-activity-weighted specificity 0-1; fa>0.3 PAMs (drives score/color)"),
    ]
    ot_descs = {
        "enpam_gb_score": "EnPAM-GB (Kim 2020) enCas12a on-target efficiency",
        "encas12a_hf1_score": "seq-DeepCpf1 enAsCas12a-HF1 variant (Chen 2025) higher-fidelity on-target efficiency",
    }
    for c in ontarget_cols:
        extra.append((c, ot_descs.get(c, "enCas12a on-target efficiency")))
    width = max(len(f) for f, _ in extra) + 2
    body = [f'    string  {(f + ";").ljust(width)}"{d}"' for f, d in extra]
    tail = [
        '    string  _mouseOver;   "Hover label (mouseOver)"',
        '    bigint  _offset;      "Uncompressed byte offset into crisprDetails.tab, or 0 if non-unique"',
        "    )",
        "",
    ]
    return "\n".join(head + body + tail)


def _write_track(
    out: dict,
    as_content: str,
    outdir: str,
    details_path: str,
    bb_name: str,
    as_name: str,
    chrom_sizes: Optional[str],
    run_tools: bool,
) -> Dict[str, str]:
    """Shared tail for :func:`build_track`: write the sorted BED + autoSql, then (if
    the UCSC tools are on PATH) bgzip+index the details and build the bigBed. ``out``
    is the ordered column dict; its key order MUST match ``as_content``."""
    bed = pd.DataFrame(out)
    bed_path = os.path.join(outdir, "guides.bed")
    bed.sort_values(["#chr", "start"]).to_csv(bed_path, sep="\t", index=False, header=False)
    as_path = os.path.join(outdir, as_name)
    with open(as_path, "w") as fh:
        fh.write(as_content)
    artifacts = {"bed": bed_path, "details_tab": details_path, "autosql": as_path}
    if run_tools:
        bgzip = _have_tool("bgzip")
        if bgzip:
            subprocess.run([bgzip, "-f", "-i", details_path], check=True)
            artifacts["details_gz"] = details_path + ".gz"
            artifacts["details_gzi"] = details_path + ".gz.gzi"
        bedtobigbed = _have_tool("bedToBigBed")
        if bedtobigbed and chrom_sizes:
            bb_path = os.path.join(outdir, bb_name)
            subprocess.run(
                [bedtobigbed, "-tab", "-type=bed9+", f"-as={as_path}", bed_path, chrom_sizes, bb_path],
                check=True,
            )
            artifacts["bigbed"] = bb_path
    return artifacts


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
    genome_fa: Optional[str] = None,
    enzyme: Optional[str] = None,
    bb_name: Optional[str] = None,
    details_name: Optional[str] = None,
    as_name: Optional[str] = None,
    encas12a_weights: Optional[str] = None,
    encas12a_list_floor: float = 0.023,
) -> Dict[str, str]:
    """Assemble the UCSC track from ``guide_df`` + crispr-ots output at
    ``enum_prefix`` (Mode-1 = ``enum_prefix``, Mode-2 = ``enum_prefix.ot.tsv``,
    sidecar = ``enum_prefix.guides.tsv``).

    ``guide_df`` must have: ``id, chrom, start, stop, strand, guideSeq, pam`` and
    any of ``DISPLAY_SCORE_COLS``. ``start``/``stop`` bound the 23-nt spacer
    (0-based, + strand). Returns a dict of the artifacts written.

    ``genome_fa`` (an indexed FASTA, ``.fai`` alongside) enables the AsCas12a
    TTTV correction: crispr-ots scores the AsCas12a (2xNLS) CFD PAM-blind, so
    without the reference PAM a TTTT off-target (which AsCas12a cannot cleave)
    would wrongly count against the AsCas12a specificity. With the genome, TTTT
    off-target CFDs are zeroed in the details popup (Layer 2) and the AsCas12a
    close-off-target grey is driven by TTTV-restricted mismatch counts (Layer 3).
    """
    # Per-enzyme filenames + on-target columns. enzyme="ascas12a" emits a single-matrix
    # AsCas12a track (TTTV-only off-target search, no TTTN/enCas12a); enzyme=None keeps
    # the legacy dual (enCas12a-primary + AsCas12a-secondary) track.
    if enzyme == "ascas12a":
        bb_name = bb_name or "ascas12a.bb"
        details_name = details_name or "ascas12a_crisprDetails.tab"
        as_name = as_name or "ascas12aTargets.as"
        score_cols = [c for c in ["deepcpf1_score", "ascas12a_deepcpf1_score"] if c in guide_df.columns]
    elif enzyme == "encas12a":
        bb_name = bb_name or "encas12a.bb"
        details_name = details_name or "encas12a_crisprDetails.tab"
        as_name = as_name or "encas12aTargets.as"
        score_cols = [c for c in ["enpam_gb_score", "encas12a_hf1_score"] if c in guide_df.columns]
    else:
        bb_name = bb_name or BB_NAME
        details_name = details_name or DETAILS_NAME
        as_name = as_name or AS_NAME
        score_cols = [c for c in DISPLAY_SCORE_COLS if c in guide_df.columns]

    os.makedirs(outdir, exist_ok=True)
    mode1 = pd.read_csv(enum_prefix)
    sidecar_path = enum_prefix + ".guides.tsv"
    mode2_path = enum_prefix + ".ot.tsv"

    # --- crisprDetails.tab + offsets ---
    # A genome enables the AsCas12a-TTTV correction (Layer 2 + Layer 3) for the legacy
    # dual track; the ascas12a path uses a TTTV-only index (no TTTT off-targets exist),
    # so no genome is needed and Mode-1 counts are already TTTV.
    genome = Genome(genome_fa) if (genome_fa and enzyme != "ascas12a") else None

    # enCas12a: fold the four guidescan specificities into the SAME Mode-2 stream. Needs the
    # reference PAM per off-target (genome) + on_count per guide (spec denominator) + the
    # PAM-activity weights. Scan UNFLOORED (--cfd-threshold 0) so the spec sums are exact;
    # encas12a_list_floor floors only the displayed popup. Specs land in enc_specs[id].
    enc_weights = enc_oncount = None
    enc_specs: Dict[str, tuple] = {}
    if enzyme == "encas12a":
        if genome is None:
            raise ValueError(
                "\n\tenCas12a track requires the reference FASTA (genome_fa, indexed) to "
                "reconstruct each off-target's PAM for the activity-weighted specificity.\n"
            )
        enc_weights = load_encas12a_weights(encas12a_weights)
        _side = pd.read_csv(sidecar_path, sep="\t")
        enc_oncount = (
            dict(zip(_side["id"], pd.to_numeric(_side["on_count"], errors="coerce").fillna(1).astype(int)))
            if "on_count" in _side.columns
            else {}
        )

    details_path = os.path.join(outdir, details_name)
    offsets, tttv_mm = _build_details(
        mode2_path,
        sidecar_path,
        mode1,
        details_path,
        list_cap,
        blank_threshold,
        genome=genome,
        weights=enc_weights,
        oncount=enc_oncount,
        specs_out=(enc_specs if enzyme == "encas12a" else None),
        list_floor=(encas12a_list_floor if enzyme == "encas12a" else 0.0),
    )

    # --- join guides with Mode-1 aggregates on the composite id ---
    # A genomic site is one track feature. Some multi-PAM guide sets (e.g. enCas12a) emit a
    # site as an EXACT-duplicate row, so both the guide table and Mode-1 can carry duplicate
    # ids (identical off-target aggregates). Keep the first of each id -- this makes the id
    # index unique for the join and draws each site once. AsCas12a/TTTV ids are already unique
    # (no-op there).
    m1 = mode1.drop_duplicates(subset="id", keep="first").set_index("id")
    g = guide_df.drop_duplicates(subset="id", keep="first").reset_index(drop=True)
    for col in ["specificity_tttv", "specificity_tttn", "off_target_count", "mismatch_counts", "dropped"]:
        g[col] = g["id"].map(m1[col]) if col in m1.columns else np.nan
    dropped = g["dropped"].fillna(0).to_numpy() == 1

    # Single-pass dual scoring: the secondary matrix's specificity rides in the
    # SAME Mode-1, as `specificity_tttn_<matrix>` / `specificity_tttv_<matrix>`.
    #
    # The suffix is derived from the matrix name (e.g. `2xnls_cas12a`,
    # `encas12a`). It used to be the hardcoded literal `2xnls` regardless of
    # which matrix was actually second, so matching on that literal would both
    # miss the current column names and, historically, mislabel the enzyme
    # whenever the two `--score` matrices were passed in the other order.
    # Detect by prefix instead, and keep accepting the old literal so Mode-1
    # files generated before the rename still load.
    extra_spec_data = []  # (field_label, display_name, spec_tttv, spec_tttn)
    ids = g["id"]
    for suffix in secondary_specificity_suffixes(m1.columns):
        label = SECONDARY_MATRIX_LABELS.get(suffix, suffix)
        tttv_col = f"specificity_tttv_{suffix}"
        tttn_col = f"specificity_tttn_{suffix}"
        etttv = (
            pd.to_numeric(ids.map(m1[tttv_col]), errors="coerce")
            if tttv_col in m1.columns
            else pd.Series(np.nan, index=g.index)
        )
        etttn = (
            pd.to_numeric(ids.map(m1[tttn_col]), errors="coerce")
            if tttn_col in m1.columns
            else pd.Series(np.nan, index=g.index)
        )
        extra_spec_data.append((label, label, etttv, etttn))

    # --- coordinates: spacer BED -> 27-nt site (5'-PAM, strand-aware) ---
    start = pd.to_numeric(g["start"], errors="coerce").fillna(0).astype(int).to_numpy()
    stop = pd.to_numeric(g["stop"], errors="coerce").fillna(0).astype(int).to_numpy()
    is_plus = g["strand"].astype(str).to_numpy() == "+"
    chrom_start = np.maximum(np.where(is_plus, start - PAM_LEN, start), 0)
    chrom_end = np.where(is_plus, stop, stop + PAM_LEN)
    thick_start = np.where(is_plus, chrom_start + PAM_LEN, chrom_start)
    thick_end = np.where(is_plus, chrom_end, chrom_end - PAM_LEN)

    # unique: non-dropped guides cleared the --threshold 0 screen (no 0-mm off-target);
    # dropped (non-unique) -> False.
    unique = np.where(dropped, "False", "True")

    # --- AsCas12a single-enzyme track (TTTV-only off-target search) ---
    # One matrix (cfd-cas12a:2xnls), one specificity (TTTV == TTTN under a TTTV index).
    # Two AsCas12a on-target models are shown: DeepCpf1 (Kim 2018) and the AsCas12a
    # seq-DeepCpf1 variant (Chen 2025). No enCas12a, no TTTN references, no genome.
    if enzyme == "ascas12a":
        spec = pd.to_numeric(g.get("specificity_tttv"), errors="coerce")
        spec_pct = _pct(spec).to_numpy()
        score = np.clip(np.rint(spec.fillna(0).to_numpy() * 100), 0, 1000).astype(int)
        # percentile per on-target model (score_cols is already the AsCas12a-only set)
        ot_pct = {c: _pct(g[c]).to_numpy() for c in score_cols}
        OT_LABELS = {"deepcpf1_score": "DeepCpf1", "ascas12a_deepcpf1_score": "AsCas12a-DeepCpf1"}
        # The color ramp is driven by ONE on-target model: DeepCpf1 (Kim 2018), else the
        # AsCas12a-variant model. Both are still emitted as display columns + hover.
        ramp_col = next((c for c in ("deepcpf1_score", "ascas12a_deepcpf1_score") if c in ot_pct), None)
        ramp_pct = ot_pct[ramp_col] if ramp_col else np.full(len(g), np.nan)
        # Mode-1 mismatch_counts are already TTTV (TTTV-only index), so the close-off-target
        # grey uses them directly -- no genome re-derivation.
        buckets = _mismatch_buckets(g["mismatch_counts"])
        mm0, mm1, mm2 = buckets["mm0"].to_numpy(), buckets["mm1"].to_numpy(), buckets["mm2"].to_numpy()
        # Three colorings, as Cas9 ships: itemRgb is the blended default (grey by
        # specificity, else the efficiency ramp on the PRIMARY model);
        # _crisprDeepCpf1Color is the same scheme on the SECOND on-target model (Cas9
        # fills this slot with Moreno-Mateos and calls it _crisprScanColor, which would
        # mislabel a Cas12a model, so the field is named for what it actually holds);
        # _specColor is specificity alone.
        item_rgb = _item_rgb(dropped, mm0, mm1, mm2, spec.to_numpy(), ramp_pct)
        alt_col = next((c for c in score_cols if c != ramp_col), None)
        alt_pct = ot_pct[alt_col] if alt_col else np.full(len(g), np.nan)
        deepcpf1_alt_color = _item_rgb(dropped, mm0, mm1, mm2, spec.to_numpy(), alt_pct)
        spec_color = _spec_color(dropped, mm0, mm1, mm2, spec.to_numpy())

        def _hover_as(i):
            if dropped[i]:
                parts = ["Not unique"]
            else:
                p = spec_pct[i]
                parts = [f"AsCas12a Spec: {int(round(p))}" if not np.isnan(p) else "AsCas12a Spec: 100"]
            for c in score_cols:
                pc = ot_pct[c][i]
                if not np.isnan(pc):
                    parts.append(f"{OT_LABELS.get(c, c)}: {int(round(pc))}")
            return ", ".join(parts)

        out = {
            "#chr": g["chrom"].to_numpy(),
            "start": chrom_start,
            "stop": chrom_end,
            "name": "",
            "score": score,
            "strand": g["strand"].to_numpy(),
            "thickStart": thick_start,
            "thickEnd": thick_end,
            "itemRgb": item_rgb,
            "_crisprDeepCpf1Color": deepcpf1_alt_color,
            "_specColor": spec_color,
            "guideSeq": g["guideSeq"].to_numpy(),
            "pam": g["pam"].to_numpy(),
            "unique": unique,
            "AsCas12a_specificity": _fmt_raw(spec).to_numpy(),
        }
        # Raw and percentile as SEPARATE numeric fields, as Cas9 does (fusi/fusiPerc).
        # Packing them into one "88% (70.25)" string reads fine on the details page but
        # leaves neither value filterable in the Table Browser or the API.
        for c in score_cols:
            base = _cas9_style_name(c)
            out[base] = _fmt_raw(pd.to_numeric(g[c], errors="coerce")).to_numpy()
            out[base + "Perc"] = [
                "" if np.isnan(v) else str(int(round(v))) for v in ot_pct[c]
            ]
        out["_mouseOver"] = [_hover_as(i) for i in range(len(g))]
        out["_offset"] = [offsets.get(gid, 0) for gid in g["id"].to_numpy()]
        as_content = build_autosql_ascas12a(score_cols)
        return _write_track(out, as_content, outdir, details_path, bb_name, as_name, chrom_sizes, run_tools)

    # --- enCas12a single-enzyme track (broad fa>0.3 off-target search) ---
    # Four guidescan specificities (raw + PAM-activity-weighted, at the fa>0.6 and fa>0.3
    # cutoffs) folded in from the Mode-2 stream (enc_specs). Column 5 / score / low-spec grey
    # use the fa>0.3 PAM-weighted specificity; the efficiency ramp uses EnPAM-GB. Two on-target
    # models: EnPAM-GB (exact) + seq-DeepCpf1 enAsCas12a-HF1 (higher-fidelity proxy). enCas12a
    # cleaves TTTT, so the close-off-target grey uses the (unrestricted) Mode-1 counts directly.
    if enzyme == "encas12a":
        n = len(g)
        # Mode-1 mismatch buckets (enCas12a cleaves TTTT, so no TTTV restriction). A guide with
        # a perfect (0-mm) off-target elsewhere is non-unique: its specificity is degenerate (a
        # repeat whose off-target sits at a non-fa0.6 PAM would even read spec_fa0.6 = 1.0), so
        # blank it like a dropped guide and mark unique=False. mm1/mm2 drive close-off grey.
        buckets = _mismatch_buckets(g["mismatch_counts"])
        mm0, mm1, mm2 = buckets["mm0"].to_numpy(), buckets["mm1"].to_numpy(), buckets["mm2"].to_numpy()
        nonuniq = dropped | (mm0 > 0)

        s06 = np.full(n, np.nan)
        s06w = np.full(n, np.nan)
        s03 = np.full(n, np.nan)
        s03w = np.full(n, np.nan)
        for i, gid in enumerate(g["id"].to_numpy()):
            if nonuniq[i]:  # non-unique: no operative specificity (greyed, blank columns)
                continue
            rec = enc_specs.get(gid)
            if rec is not None:
                s06[i], s06w[i], s03[i], s03w[i] = rec[0], rec[1], rec[2], rec[3]
            else:  # unique guide with no off-target -> 1/on_count (1.0 when on_count == 1)
                v = 1.0 / max(1, (enc_oncount or {}).get(gid, 1))
                s06[i] = s06w[i] = s03[i] = s03w[i] = v
        score_spec = s03w  # fa>0.3 PAM-weighted specificity drives col5/score + low-spec grey
        score = np.clip(np.rint(np.nan_to_num(score_spec, nan=0.0) * 100), 0, 1000).astype(int)
        unique_en = np.where(nonuniq, "False", "True")

        enpam_pct = _pct(g["enpam_gb_score"]).to_numpy() if "enpam_gb_score" in g.columns else np.full(n, np.nan)
        item_rgb = _item_rgb(dropped, mm0, mm1, mm2, score_spec, enpam_pct)

        ot_cols = [c for c in ("enpam_gb_score", "encas12a_hf1_score") if c in g.columns]
        OT_LABELS = {"enpam_gb_score": "EnPAM-GB", "encas12a_hf1_score": "enAsCas12a-HF1"}
        ot_pct = {c: _pct(g[c]).to_numpy() for c in ot_cols}

        def _hover_en(i):
            if nonuniq[i]:
                parts = ["Not unique"]
            else:
                sv = score_spec[i]
                parts = ["enCas12a Spec (fa>0.3 wt): " + (f"{sv:.3f}" if not np.isnan(sv) else "1.000")]
            for c in ot_cols:
                pc = ot_pct[c][i]
                if not np.isnan(pc):
                    parts.append(f"{OT_LABELS.get(c, c)}: {int(round(pc))}")
            return ", ".join(parts)

        out = {
            "#chr": g["chrom"].to_numpy(),
            "start": chrom_start,
            "stop": chrom_end,
            "name": "",
            "score": score,
            "strand": g["strand"].to_numpy(),
            "thickStart": thick_start,
            "thickEnd": thick_end,
            "itemRgb": item_rgb,
            "guideSeq": g["guideSeq"].to_numpy(),
            "pam": g["pam"].to_numpy(),
            "unique": unique_en,
            "enCas12a_spec_fa06": _fmt_raw(pd.Series(s06)).to_numpy(),
            "enCas12a_spec_fa06_wt": _fmt_raw(pd.Series(s06w)).to_numpy(),
            "enCas12a_spec_fa03": _fmt_raw(pd.Series(s03)).to_numpy(),
            "enCas12a_spec_fa03_wt": _fmt_raw(pd.Series(s03w)).to_numpy(),
        }
        for c in ot_cols:
            out[c] = _fmt_pct_raw(g[c], ot_pct[c]).to_numpy()
        out["_mouseOver"] = [_hover_en(i) for i in range(n)]
        out["_offset"] = [offsets.get(gid, 0) for gid in g["id"].to_numpy()]
        as_content = build_autosql_encas12a(ot_cols)
        return _write_track(out, as_content, outdir, details_path, bb_name, as_name, chrom_sizes, run_tools)

    # --- scores / colors / display strings ---
    # enCas12a display specificity (primary matrix). enCas12a DOES cleave TTTT, so its
    # operative aggregate is TTTN; both TTTV and TTTN are emitted as display columns.
    spec_tttv = pd.to_numeric(g.get("specificity_tttv"), errors="coerce")
    spec_tttn = pd.to_numeric(g.get("specificity_tttn"), errors="coerce")
    spec_pct = _pct(spec_tttn).to_numpy()
    # Column 5 (score), the UCSC score filter, and the low-spec greying use AsCas12a (WT)
    # **TTTV** specificity when the 2xNLS matrix is present, else EnCas12a TTTN. AsCas12a
    # does not cleave TTTT, so its operative specificity excludes TTTT off-targets.
    # extra_spec_data rows are (field, display, spec_tttv, spec_tttn) -> [0][2] = AsCas12a TTTV.
    score_spec = extra_spec_data[0][2] if extra_spec_data else spec_tttv
    score = np.clip(np.rint(score_spec.fillna(0).to_numpy() * 100), 0, 1000).astype(int)
    # Greying keys off the SAME raw AsCas12a-TTTV specificity as the score, via an ABSOLUTE
    # two-tier threshold (not a percentile). Hover still reports EnCas12a Spec via spec_pct.
    color_spec_raw = score_spec.to_numpy()

    # Close-off-target grey (120,120,120) is driven by 1/2-mm counts. With a genome
    # these are TTTV-restricted (Layer 3): AsCas12a does not cleave TTTT, so a close
    # TTTT off-target must not grey the guide. Mode-1 mismatch_counts are TTTN (they
    # include TTTT); without a genome, fall back to them.
    if tttv_mm:
        tv = np.array([tttv_mm.get(gid, (0, 0, 0)) for gid in g["id"].to_numpy()], dtype=int)
        mm0, mm1, mm2 = tv[:, 0], tv[:, 1], tv[:, 2]
    else:
        buckets = _mismatch_buckets(g["mismatch_counts"])
        mm0 = buckets["mm0"].to_numpy()
        mm1 = buckets["mm1"].to_numpy()
        mm2 = buckets["mm2"].to_numpy()
    enpam_pct = _pct(g["enpam_gb_score"]).to_numpy() if "enpam_gb_score" in g.columns else np.full(len(g), np.nan)
    enseq_pct = (
        _pct(g["enseq_deepcpf1_score"]).to_numpy() if "enseq_deepcpf1_score" in g.columns else np.full(len(g), np.nan)
    )
    has_ascas12a = "ascas12a_deepcpf1_score" in g.columns
    ascas12a_pct = _pct(g["ascas12a_deepcpf1_score"]).to_numpy() if has_ascas12a else np.full(len(g), np.nan)
    # itemRgb efficiency ramp is driven by AsCas12a-DeepCpf1 (falls back to
    # EnCas12a-DeepCpf1 when the AsCas12a score isn't present). Also shown on hover.
    color_pct = ascas12a_pct if has_ascas12a else enseq_pct
    item_rgb = _item_rgb(dropped, mm0, mm1, mm2, color_spec_raw, color_pct)

    out = {
        "#chr": g["chrom"].to_numpy(),
        "start": chrom_start,
        "stop": chrom_end,
        "name": "",
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

    # each extra matrix's OPERATIVE specificity percentile for the hover: AsCas12a is
    # TTTV (excludes TTTT off-targets, which it does not cleave).
    extra_spec_pct = [(display, _pct(etttv).to_numpy()) for _field, display, etttv, _etttn in extra_spec_data]

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

    # --- write BED (sorted), autoSql, tools (shared tail) ---
    spec_matrices = [("enCas12a", "EnCas12a")] + [(field, display) for field, display, _tv, _tn in extra_spec_data]
    as_content = build_autosql(score_cols, spec_matrices)
    return _write_track(out, as_content, outdir, details_path, bb_name, as_name, chrom_sizes, run_tools)


# ----------------------------------------------------------------------------
# genome-wide percentile re-grading (chunked builds only)
# ----------------------------------------------------------------------------
# A chunked build scores each ~chunk_max slice of the guide BED as an isolated
# run, so every percentile it computes is a rank WITHIN that slice. The guide BED
# is position-ordered, so a slice is a contiguous stretch of one or two
# chromosomes: a percentile would describe a guide against its regional
# neighbours rather than against the genome, and the efficiency ramp colours
# derived from one would mean different things in different places. That is wrong
# for a genome-wide track, so the merge re-ranks every percentile across all
# chunks before building the bigBed.
#
# This is possible because the raw score columns are exact per guide: the ranks
# are recoverable from them with no rescan and no re-scoring. Only ranked values
# need it. ``score`` (specificity x100) and ``_specColor`` (absolute cutoffs on
# the raw specificity) are already genome-wide and are left alone.
# bytes, not str: the re-grade works on raw BED lines, and a str set here would
# never match, silently replacing every grey with a ramp colour.
_RAMP_GREYS = frozenset({b"150,150,150", b"120,120,120", b"90,90,90", b"50,50,50"})

# Per-enzyme map of what to re-grade. ``models`` is one (raw field, hover label,
# colour field) triple per on-target model, in the order build_track writes the
# hover; the first model's colour is the default ``itemRgb``.
REGRADE_SPECS: Dict[str, dict] = {
    "ascas12a": {
        "unique": "unique",
        "spec": "AsCas12a_specificity",
        "spec_label": "AsCas12a Spec",
        "hover": "_mouseOver",
        "models": [
            ("deepcpf1", "DeepCpf1", "itemRgb"),
            ("ascas12aDeepcpf1", "AsCas12a-DeepCpf1", "_crisprDeepCpf1Color"),
        ],
    },
}


def autosql_field_names(as_path: str) -> list[str]:
    """Ordered field names declared by an autoSql file, i.e. the BED column order."""
    import re

    names = []
    with open(as_path) as fh:
        for line in fh:
            m = re.match(r"\s+(?:string|uint|char|bigint)\S*\s+(\w+);", line)
            if m:
                names.append(m.group(1))
    return names


def _rank_pct(counts: Dict[bytes, int]) -> Dict[bytes, float]:
    """``{raw value -> percentile}`` from a value histogram, reproducing
    ``_pct`` (pandas ``rank(pct=True, method="average") * 100``, 2 dp) over the
    non-blank values. Ranking the histogram instead of the values keeps this
    O(distinct) in memory: the raw columns are written to 4 decimals, so a
    genome of 10^8 guides still has at most ~10^6 distinct keys."""
    n = sum(counts.values())
    out: Dict[bytes, float] = {}
    if not n:
        return out
    below = 0
    for v in sorted(counts, key=float):
        c = counts[v]
        out[v] = round(100.0 * (below + (c + 1) / 2.0) / n, 2)
        below += c
    return out


def _ramp_colour(pct: Optional[float]) -> bytes:
    """The efficiency ramp. A missing percentile is treated as fully efficient,
    matching :func:`_spec_eff_color`."""
    enp = 100.0 if pct is None else pct
    if enp <= EFF_RAMP_LOW:
        return b"255,100,100"
    return b"255,255,0" if enp <= EFF_RAMP_HIGH else b"0,200,0"


def _regrade_bed_percentiles(
    bed_in: str,
    bed_out: str,
    fields: Sequence[str],
    spec: dict,
    counts: Dict[str, Dict[bytes, int]],
) -> int:
    """Rewrite ``bed_in`` to ``bed_out`` with every percentile re-ranked over the
    whole genome: the ``...Perc`` columns, the hover percentiles, and the
    efficiency-ramp colours.

    A colour cell that is already a grey is left untouched. The greys encode
    specificity and uniqueness, which :func:`_item_rgb` applies after the ramp
    and which are genome-wide to begin with, so a non-grey cell is always a pure
    ramp colour and is the only thing to re-derive."""
    idx = {f: i for i, f in enumerate(fields)}
    lut = {f: _rank_pct(c) for f, c in counts.items()}
    spec_lut = lut.get(spec["spec"], {})
    i_spec = idx[spec["spec"]]
    i_uniq = idx[spec["unique"]]
    i_hover = idx[spec["hover"]]
    label_spec = spec["spec_label"].encode()
    models = [
        (idx[raw], idx.get(raw + "Perc"), idx[col], label.encode(), lut.get(raw, {}))
        for raw, label, col in spec["models"]
    ]

    n = 0
    with open(bed_in, "rb") as fin, open(bed_out, "wb") as fout:
        for row in fin:
            row = row.rstrip(b"\n")
            if not row:
                continue
            f = row.split(b"\t")
            # hover: "Not unique" for a non-unique guide, else the specificity
            # percentile, then one entry per model that scored.
            if f[i_uniq] == b"False":
                parts = [b"Not unique"]
            else:
                sp = spec_lut.get(f[i_spec])
                parts = [label_spec + b": " + b"%d" % round(100.0 if sp is None else sp)]
            for i_raw, i_pct, i_col, label, model_lut in models:
                pct = model_lut.get(f[i_raw]) if f[i_raw] else None
                if i_pct is not None:
                    f[i_pct] = b"" if pct is None else b"%d" % round(pct)
                if f[i_col] not in _RAMP_GREYS:
                    f[i_col] = _ramp_colour(pct)
                if f[i_raw]:
                    parts.append(label + b": " + b"%d" % round(pct))
            f[i_hover] = b", ".join(parts)
            fout.write(b"\t".join(f) + b"\n")
            n += 1
    return n


def merge_ucscgb_tracks(
    track_dirs: Sequence[str],
    out: str,
    chrom_sizes: str,
    as_file: str,
    bedtobigbed: str = "bedToBigBed",
    sort_mem: str = "8G",
    sort_threads: str = "8",
    keep_bed: bool = False,
    details_name: str = DETAILS_NAME,
    bb_name: str = BB_NAME,
    as_name: str = AS_NAME,
    enzyme: str = "ascas12a",
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

    Every percentile is then re-ranked across all chunks (see
    :func:`_regrade_bed_percentiles`), since each chunk could only rank within
    itself. ``enzyme`` selects which columns that applies to; an enzyme with no
    entry in :data:`REGRADE_SPECS` is merged without re-grading and warns, because
    its percentiles would stay chunk-local.

    Returns the artifact paths.
    """
    header = b"_mismatchCounts\t_crisprOfftargets\n"

    # Columns whose percentiles the merge must re-rank genome-wide, and the raw
    # columns to histogram while streaming the per-chunk BEDs.
    spec = REGRADE_SPECS.get(enzyme)
    as_fields = autosql_field_names(as_file)
    if spec is not None:
        needed = [spec["spec"], spec["unique"], spec["hover"]]
        needed += [raw for raw, _, _ in spec["models"]] + [col for _, _, col in spec["models"]]
        missing = [f for f in needed if f not in as_fields]
        if missing:
            raise ValueError(f"{as_file} is missing fields needed to re-grade percentiles: {missing}")
    else:
        print(
            f"\n\tWARNING: no percentile re-grading rule for enzyme={enzyme!r}; the merged track's "
            f"percentiles and efficiency-ramp colours stay ranked WITHIN each chunk, not genome-wide.\n",
            flush=True,
        )
    count_fields = [spec["spec"]] + [raw for raw, _, _ in spec["models"]] if spec else []
    count_idx = {f: as_fields.index(f) for f in count_fields}
    counts: Dict[str, Dict[bytes, int]] = {f: {} for f in count_fields}

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
    merged_tab = os.path.join(out, details_name)
    merged_bed = os.path.join(out, "guides.merged.bed")

    n_guides = n_unique = 0
    running = 0  # byte offset in the merged uncompressed details
    with open(merged_tab, "wb") as outf, open(merged_bed, "wb") as bedout:
        outf.write(header)
        running = len(header)
        for tdir in track_dirs:
            det = os.path.join(tdir, details_name + ".gz")
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
                    # histogram the raw columns while the row is already split, so
                    # the genome-wide ranks cost no extra pass over the chunk BEDs
                    for f_name, f_i in count_idx.items():
                        v = fields[f_i]
                        if v:
                            counts[f_name][v] = counts[f_name].get(v, 0) + 1
                    bedout.write(b"\t".join(fields) + b"\n")
                    n_guides += 1
            offmap.clear()
            print(f"  merged {os.path.basename(tdir)}: cum_guides={n_guides} merged_bytes={running}", flush=True)

    # 3. re-rank every percentile across all chunks (each chunk could only rank
    # within itself), then sort the merged BED and build the bigBed
    if spec is not None:
        regraded_bed = merged_bed + ".regraded"
        n_re = _regrade_bed_percentiles(merged_bed, regraded_bed, as_fields, spec, counts)
        if n_re != n_guides:
            raise RuntimeError(f"percentile re-grade wrote {n_re} rows, expected {n_guides}")
        os.replace(regraded_bed, merged_bed)
        print(
            "  re-graded percentiles genome-wide over "
            + ", ".join(f"{f} ({len(counts[f])} distinct)" for f in count_fields),
            flush=True,
        )

    sorted_bed = os.path.join(out, "guides.sorted.bed")
    subprocess.run(
        f"LC_ALL=C sort -S{sort_mem} --parallel={sort_threads} -T {out} -k1,1 -k2,2n {merged_bed} > {sorted_bed}",
        shell=True,
        check=True,
    )
    bb = os.path.join(out, bb_name)
    subprocess.run([bedtobigbed, "-tab", "-type=bed9+", f"-as={as_file}", sorted_bed, chrom_sizes, bb], check=True)
    # 4. bgzip + index the merged details (offsets are uncompressed; .gzi handles the seek)
    subprocess.run(["bgzip", "-f", "-@", "8", merged_tab], check=True)
    subprocess.run(["bgzip", "-r", merged_tab + ".gz"], check=True)
    shutil.copy(as_file, os.path.join(out, as_name))
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
        "autosql": os.path.join(out, as_name),
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
