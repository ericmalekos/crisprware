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


# ----------------------------------------------------------------------------
# autoSql
# ----------------------------------------------------------------------------
def _pct(series: pd.Series) -> pd.Series:
    """Percentile rank in [0, 100] (higher value -> higher percentile)."""
    x = pd.to_numeric(series, errors="coerce")
    return (x.rank(pct=True, method="average") * 100).round(2)


# Displayed precision of the specificity columns. 3 places, not 4: a CFD aggregate
# does not support the fourth, and the Cas9 track shows none. This is safe to change
# where the on-target precision is not, because specificity is never ranked, so no
# percentile moves with it. At 3 places the displayed value differs from the
# full-precision float by at most 0.0005, i.e. 0.05 of a `score` unit.
SPEC_DECIMALS = 3


def _pct_of_displayed(raw: pd.Series, decimals: int) -> pd.Series:
    """Percentile rank of a cleavage score AS THE TRACK DISPLAYS IT, i.e. of the
    rounded value :func:`_fmt_score` writes, not of the full-precision float behind it.

    This is what makes a percentile reproducible from the published file, and it
    is what lets a chunked build re-rank percentiles genome-wide at merge time
    (:func:`_regrade_bed_percentiles`), which only has the written columns to work
    from. Ranking the hidden precision instead silently makes the two build paths
    disagree: on sacCer3, 711 guides whose specificity rounds to a displayed 1.0
    rank near the 5th percentile on full floats but join the 53rd-percentile tie
    once rounded, so 4114 of 590,425 hover labels differed between an unchunked
    build and a chunked one. Ranking the displayed value removes that whole class
    of disagreement rather than narrowing it.

    It follows that coarsening the display RE-RANKS, so reducing a model's decimals
    is never a cosmetic change. Ranking the 4-decimal value while showing an integer
    would leave 68 of 112 displayed DeepCpf1 values carrying two or more different
    percentiles (measured over 2,120,555 rn7 chr20 guides), which is precisely what
    :func:`validate_track_bed` rejects. Re-ranking preserves the invariant; the price
    is that the displayed percentile moves for about 37% of ramp-coloured guides, by
    at most 3.31 points, and the ramp colour itself for 0.46% of them."""
    return _pct(pd.to_numeric(_fmt_score(raw, decimals), errors="coerce"))


def _fmt_raw(raw: pd.Series, decimals: int = SPEC_DECIMALS) -> pd.Series:
    """Specificity to <=``decimals`` places as a string; blank for NaN. Used ONLY for
    the specificity columns: a percentile rank is degenerate there because the vast
    majority of guides tie at a perfectly-specific 1.0, so the raw CFD specificity
    (0-1) is what's shown, and it is never ranked.

    Variable width (``round()``, so ``1.0`` not ``1.000`` and ``0.16`` not ``0.160``):
    trailing zeros are stripped but one is always kept, so an exact 1 stays ``1.0``.

    Unlike :func:`_fmt_score`, rounding here is genuinely display-only, because no
    percentile is derived from specificity. It does narrow the agreement between
    ``score`` and the displayed value, but only to 0.05 score units at 3 places,
    well inside the tolerance :func:`validate_track_bed` allows."""
    out = []
    for r in raw:
        out.append("" if pd.isna(r) else f"{round(float(r), decimals)}")
    return pd.Series(out, index=raw.index, dtype="object")


def _fmt_score(raw: pd.Series, decimals: int) -> pd.Series:
    """An on-target cleavage score as the track displays it: fixed-width to
    ``decimals`` places, blank for NaN.

    The precision is per model and is declared in :data:`ENZYME_TRACK_SPECS`, because
    it is a property of the model's output scale, not of the enzyme. ``50.9896`` is
    more precision than any of these models supports and more than the Cas9 track
    shows, which reports all six of its score fields as plain integers. Rounding puts
    each model in the same granularity band as Cas9's 92 to 126 distinct values:
    DeepCpf1, which runs 0.96 to 114.7, gives 112 distinct integers, while the three
    0-to-1 models (seq-DeepCpf1, EnPAM-GB, enAsCas12a-HF1) give 92 to 94 distinct
    values at 2 decimals and a useless 2 as integers.

    Fixed width rather than ``round()``'s variable width, so the column reads
    consistently: ``0.70`` beside ``0.73``, not ``0.7``."""
    out = []
    for r in raw:
        out.append("" if pd.isna(r) else f"{float(r):.{decimals}f}")
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
SPEC_GREY_BAD = 0.20  # raw specificity < 0.20       -> grey (bad)
# There is deliberately no middle tier. A 0.20-to-0.50 guide used to be greyed
# "90,90,90", which hid its efficiency ramp for a specificity that is poor but still
# usable; it now shows the ramp like any other specific guide and its specificity is
# read from `score`, `_specColor` and the hover instead.

# Efficiency-percentile cuts for the red/yellow/green ramp, taken from the Cas9
# crisprAllTargets track description so the two tracks read the same way:
# "low predicted cleavage: ... percentile <= 30", "medium ...: > 30 and < 55",
# "high ...: > 55". (These were 50/75 before, markedly stricter than Cas9, so a
# guide at percentile 60 showed yellow here and green on the Cas9 track.)
EFF_RAMP_LOW = 30  # percentile <= 30            -> red
EFF_RAMP_HIGH = 55  # 30 < percentile <= 55      -> yellow; > 55 -> green

# Cuts for the specificity-only _specColor, on a 0-100 scale, so the raw 0-1 AsCas12a
# specificity is multiplied by 100 first.
#
# The red cut is 20, matching SPEC_GREY_BAD, and NOT the 50 read off the Cas9
# crisprAllTargets sample rows (specificity 27/37/43 -> red, 51/52/65/69 -> olive,
# 72/76 -> green). Carrying both numbers meant the track drew two different
# "bad specificity" lines on the same quantity: a guide at 0.35 kept its efficiency
# ramp in the default view, which says it is not hard to target, and read red in the
# specificity view, which says it is. One of the two had to move, and 0.20 is the one
# the default scheme, the score filter and the description page already use.
#
# This is a deliberate divergence from the Cas9 track, which greys everything at or
# below MIT 50 and so has no such split to resolve. The olive band is correspondingly
# wide (0.20 to 0.70); green above 0.70 is unchanged.
SPEC_COLOR_BAD = 20  # <= 20      -> red
SPEC_COLOR_MID = 70  # 20 to 70   -> olive; > 70 -> green

# Every track ships the same three colourings, in this order, directly after bed9.
# ``itemRgb`` is the blended default; the two alternatives are selectable from the
# browser's "Color by:" dropdown via the trackDb ``colorFields`` setting.
COLOR_FIELDS = ("itemRgb", "_crisprDeepCpf1Color", "_specColor")

# The ramp-driven colour fields, in on-target model order: the first scored model
# drives ``itemRgb``, the second drives ``_crisprDeepCpf1Color``. ``_specColor`` is
# not here because it is a function of specificity alone and is never ranked.
EFF_COLOR_FIELDS = ("itemRgb", "_crisprDeepCpf1Color")

# Per-enzyme on-target ("cleavage") models, in display order: the only thing that
# legitimately differs between enzymes. Everything else about a track -- the three
# colourings, the split raw/``...Perc`` columns, the hover layout, and the
# genome-wide percentile basis -- is universal and derived from this table, so a new
# enzyme or an extra cleavage score needs an entry here and nothing else.
#
# Each model is ``(column, hover label, displayed decimals)``. The decimals are a
# property of the model's output scale: DeepCpf1 spans roughly 1 to 115 and reads as
# an integer, while the three 0-to-1 models need 2 places (see :func:`_fmt_score` for
# the measured distinct-value counts). Changing a model's decimals re-ranks its
# percentile and can move its ramp colour, so it is a track change, not a display
# tweak.
# Hover labels follow the Cas9 track, whose mouseOver reads
# "MIT Spec. Score: 48, Doench 2016: 18%, Moreno-Mateos: 72%": the specificity label
# ends in "Spec. Score" and each cleavage figure is a percentile with a % suffix.
#
# The first AsCas12a model is cited rather than named because the obvious name is
# wrong: crisprware loads weights/Seq_deepCpf1_weights.h5, the SEQUENCE-ONLY model
# from Kim et al. 2018. "DeepCpf1" proper is the chromatin-accessibility variant and
# is not what runs, so the tracks shipped before 2026-10-19 mislabelled it.
ENZYME_TRACK_SPECS: Dict[str, dict] = {
    "ascas12a": {
        "spec_label": "AsCas12a Spec. Score",
        "score_spec": "AsCas12a_specificity",
        "ontarget": [
            ("deepcpf1_score", "Kim 2018", 0),
            ("ascas12a_deepcpf1_score", "Chen 2025", 2),
        ],
    },
    "encas12a": {
        # The qualifier stays: four specificities are reported and this names the one
        # driving score and the greys.
        "spec_label": "enCas12a Spec. Score (fa>0.3 wt)",
        "score_spec": "enCas12a_spec_fa03_wt",
        "ontarget": [
            ("enpam_gb_score", "EnPAM-GB", 2),
            ("encas12a_hf1_score", "enAsCas12a-HF1", 2),
        ],
    },
}


def _cas9_style_name(col: str) -> str:
    """``deepcpf1_score`` -> ``deepcpf1``; ``ascas12a_deepcpf1_score`` ->
    ``ascas12aDeepcpf1``. Cas9 names each efficiency model with a bare camelCase
    field plus a ``...Perc`` sibling (``fusi``/``fusiPerc``,
    ``crisprScan``/``crisprScanPerc``); this derives that base name."""
    s = col[:-6] if col.endswith("_score") else col
    head, *rest = s.split("_")
    return head + "".join(p.capitalize() for p in rest)


def _spec_eff_color(spec_raw, eff_pct) -> np.ndarray:
    """Color for a unique guide with no close (1/2-mm) off-targets: a single
    ABSOLUTE low-specificity grey below :data:`SPEC_GREY_BAD` on the raw specificity,
    else the red/yellow/green efficiency ramp on the on-target percentile. Missing
    specificity is treated as fully specific (not greyed)."""
    enp = np.where(np.isnan(eff_pct), 100.0, eff_pct)
    out = np.where(
        enp <= EFF_RAMP_LOW, "255,100,100", np.where(enp <= EFF_RAMP_HIGH, "255,255,0", "0,200,0")
    )
    sr = np.where(np.isnan(spec_raw), 1.0, spec_raw)
    out = np.where(sr < SPEC_GREY_BAD, "50,50,50", out)  # bad < 0.20
    return out


def _item_rgb(dropped, mm0, mm1, mm2, spec_raw, eff_pct) -> np.ndarray:
    """Color hierarchy (highest priority applied last): efficiency ramp / two-tier
    absolute low-specificity grey (see _spec_eff_color), then grey for close (1/2-mm)
    off-targets, then light grey for non-unique (0-mm / dropped). ``spec_raw`` is the
    raw AsCas12a-TTTV specificity. ``mm*`` are the scanner's mismatch buckets, which are
    TTTV-restricted for AsCas12a only because the index itself enumerates TTTV: the
    scanner counts every site it is given. `crot_index.assert_enzyme_pams` enforces that
    at the start of the run, because nothing downstream can tell a TTTT site apart.

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
    with green, and keep red even where ``reserved`` has already greyed out. The
    no-greys part is kept, so a greyed guide still gets a colour here, but the cuts
    are ``<= 20`` red, ``<= 70`` olive, ``> 70`` green: the red cut follows
    :data:`SPEC_GREY_BAD` rather than Cas9's 50, so the track states one
    "bad specificity" line instead of two. ``spec_raw`` is the raw AsCas12a-TTTV
    specificity (0-1).
    A non-unique guide has no meaningful specificity, so it takes the worst tier;
    a missing specificity keeps the usual "treat as fully specific" convention."""
    pct = np.where(np.isnan(spec_raw), 100.0, np.asarray(spec_raw, dtype=float) * 100.0)
    pct = np.where(dropped | (mm0 > 0), 0.0, pct)  # non-unique -> worst tier
    return np.where(
        pct <= SPEC_COLOR_BAD, "255,0,0", np.where(pct <= SPEC_COLOR_MID, "128,128,0", "0,255,0")
    )


def ordered_ontarget_cols(enzyme: str, columns: Sequence[str]) -> list[str]:
    """The enzyme's on-target columns present in ``columns``, in display order.

    Display order is what binds a model to a colour field, so it is taken from
    :data:`ENZYME_TRACK_SPECS` rather than from whatever order the scorer emitted."""
    espec = ENZYME_TRACK_SPECS.get(enzyme)
    if espec is None:
        return [c for c in columns if c.endswith("_score")]
    return [c for c, _lab, _dec in espec["ontarget"] if c in columns]


def _ontarget_block(
    g: pd.DataFrame,
    enzyme: str,
    ot_cols: Sequence[str],
    dropped,
    mm0,
    mm1,
    mm2,
    spec_raw: np.ndarray,
    score: np.ndarray,
) -> tuple:
    """The part of a track that is identical for every enzyme.

    Given the enzyme's scored on-target ("cleavage") models and the specificity that
    already drives column 5, this produces the three colourings, the split
    raw/``...Perc`` display columns and the hover label. Returns
    ``(colours, display, hover)``.

    Universal by construction, so adding a cleavage score or a new enzyme changes
    only :data:`ENZYME_TRACK_SPECS`:

    * percentiles rank the **displayed** value at each model's declared precision
      (:func:`_pct_of_displayed`), which is what lets
      :func:`_regrade_bed_percentiles` re-rank them genome-wide from the written file
      after a chunked build;
    * each model emits a raw column plus a ``...Perc`` sibling, as the Cas9 track
      does with ``fusi``/``fusiPerc``, so both halves stay filterable in the Table
      Browser instead of being packed into one ``"88% (70.25)"`` string;
    * the first model drives ``itemRgb``, the second ``_crisprDeepCpf1Color``, and
      ``_specColor`` is specificity alone;
    * the hover reports specificity as the RAW value x100, equal to ``score`` and
      consistent with ``_specColor``, never as a percentile (ranking specificity is
      degenerate: the distribution is dominated by the tie at a perfectly specific
      1.0, so a guide at 0.9997 would read as terrible while ``score`` says 100).
    """
    n = len(g)
    label = ENZYME_TRACK_SPECS.get(enzyme, {}).get("spec_label", "Spec")
    ontarget = ENZYME_TRACK_SPECS.get(enzyme, {}).get("ontarget", [])
    labels = {c: lab for c, lab, _dec in ontarget}
    decimals = {c: dec for c, _lab, dec in ontarget}

    ot_pct = {c: _pct_of_displayed(g[c], decimals[c]).to_numpy() for c in ot_cols}
    ramp = [ot_pct[c] for c in ot_cols] + [np.full(n, np.nan)] * len(EFF_COLOR_FIELDS)

    colours = {
        f: _item_rgb(dropped, mm0, mm1, mm2, spec_raw, ramp[i])
        for i, f in enumerate(EFF_COLOR_FIELDS)
    }
    colours["_specColor"] = _spec_color(dropped, mm0, mm1, mm2, spec_raw)

    display: Dict[str, np.ndarray] = {}
    for c in ot_cols:
        base = _cas9_style_name(c)
        display[base] = _fmt_score(
            pd.to_numeric(g[c], errors="coerce"), decimals[c]
        ).to_numpy()
        display[base + "Perc"] = np.array(
            ["" if np.isnan(v) else str(int(round(v))) for v in ot_pct[c]], dtype=object
        )

    nonuniq = np.asarray(dropped, dtype=bool)
    hover = []
    for i in range(n):
        parts = ["Not unique"] if nonuniq[i] else [f"{label}: {score[i]}"]
        for c in ot_cols:
            pc = ot_pct[c][i]
            if not np.isnan(pc):
                # "%" because the figure IS a percentile, as on the Cas9 track
                parts.append(f"{labels.get(c, c)}: {int(round(pc))}%")
        hover.append(", ".join(parts))

    return colours, display, hover


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
    """Write crisprDetails.tab and return ``offsets``.

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

    ``offsets`` maps ``id -> uncompressed byte offset``.

    ``genome`` is used only to read each off-target's reference 4-nt PAM for the enCas12a
    specificity weights; the scanner itself is PAM-blind. The AsCas12a TTTT exclusion that
    used to live here (zeroing a second CFD column, and a TTTV-restricted mismatch count)
    belonged to the two-matrix track and is gone with it: a single-enzyme AsCas12a track is
    scanned against a TTTV index, where no TTTT off-target exists to exclude (enforced by
    `crot_index.assert_enzyme_pams`, since a TTTN index silently reinstates every TTTT
    off-target in the counts, the popup list and the 0-mm screen), and enCas12a does
    cleave TTTT so its counts must stay unrestricted.
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
                    # Reference 4-nt PAM of this off-target (the scanner is PAM-blind),
                    # for the enCas12a specificity weights.
                    pam_str = genome.pam(p[1], pos, strand) if genome is not None else ""
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
                    # A single-enzyme track scans ONE CFD matrix, so each off-target is
                    # chrom;posStrand;scoreInt.
                    entry = f"{p[1]};{pos}{strand};{si}"
                    # list_floor floors ONLY the displayed popup (enCas12a scans unfloored so
                    # the spec sums above are exact); default 0 keeps AsCas12a/legacy behavior.
                    if cfd >= list_floor:
                        entries.append((entry, si, mm))
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
    return offsets


# ----------------------------------------------------------------------------
# secondary-matrix specificity columns
# ----------------------------------------------------------------------------
# Display names for known matrix-name slugs. crispr-ots derives the column
# suffix from the matrix's own name, so `2xNLS-Cas12a` -> `2xnls_cas12a`.
def build_autosql_ascas12a(ontarget_cols: Sequence[str] = ()) -> str:
    """autoSql for the single-enzyme **AsCas12a** track: bed9 + ``_specColor`` +
    guideSeq/pam + a single ``unique`` flag + one ``AsCas12a_specificity`` column +
    one column per on-target model in ``ontarget_cols`` + the ``_mouseOver``/``_offset``
    trailer. AsCas12a is scored against a TTTV off-target index, so there is NO
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
        '    uint    itemRgb;      "Display color R,G,B (specificity grey, else Kim et al. 2018 efficiency ramp)"',
        '    string  _crisprDeepCpf1Color; "Alternative color R,G,B by the second on-target model (Chen et al. 2025)"',
        '    string  _specColor;   "Alternative color R,G,B by AsCas12a-TTTV specificity alone"',
        '    string  guideSeq;     "Guide sequence"',
        '    string  pam;          "Protospacer Adjacent Motif (PAM)"',
    ]
    extra = [
        ("unique", "Guide has a unique on-target (TTTV off-target search)"),
        ("AsCas12a_specificity", "AsCas12a aggregate CFD specificity"),
    ]
    # Cas9's crispr.as convention exactly, as used for fusi/fusiPerc and
    # crisprScan/crisprScanPerc: "Efficiency: <citation> Score (raw|percentile)".
    # Short enough to read as a field label on its own, so it needs no "label|note"
    # split. The model is CITED rather than named because the obvious name is wrong:
    # what runs is Seq-deepCpf1, the sequence-only model, not DeepCpf1 proper.
    # No trailing whitespace in any comment: hand patches once padded them to hold a
    # byte length, that is gone, and it stays gone.
    ot_descs = {
        "deepcpf1_score": "Kim et al. 2018 Score",
        "ascas12a_deepcpf1_score": "Chen et al. 2025 Score",
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
    """autoSql for the single-enzyme **enCas12a** track: bed9 + the two alternative colour
    fields + guideSeq/pam + a unique flag + the four guidescan specificities (raw +
    PAM-activity-weighted, at the fa>0.6 and fa>0.3 activity cutoffs) + a raw and a
    ``...Perc`` column per on-target model in ``ontarget_cols`` + the ``_mouseOver``/
    ``_offset`` trailer. enCas12a is scored against the broad fa>0.3 (30-PAM) off-target
    index; the fa>0.3 PAM-weighted specificity drives the score and the greys.

    Structurally identical to :func:`build_autosql_ascas12a`: same three colourings, same
    split raw/percentile layout, same trailer. The enzymes differ only in which specificity
    columns and which cleavage models they carry. Column order MUST match the BED written by
    :func:`build_track` with ``enzyme="encas12a"``."""
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
        '    uint    itemRgb;      "Display color R,G,B (specificity grey, else EnPAM-GB efficiency ramp)"',
        '    string  _crisprDeepCpf1Color; "Alternative color R,G,B by the second on-target model (enAsCas12a-HF1, Chen et al. 2025)"',
        '    string  _specColor;   "Alternative color R,G,B by enCas12a fa>0.3 PAM-weighted specificity alone"',
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
    # Same Cas9 crispr.as convention as the AsCas12a builder; see the note there.
    ot_descs = {
        "enpam_gb_score": "EnPAM-GB, Kim et al. 2020 Score",
        "encas12a_hf1_score": "enAsCas12a-HF1, Chen et al. 2025 Score",
    }
    for c in ontarget_cols:
        d = ot_descs.get(c, "enCas12a on-target efficiency")
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


def _write_track(
    out: dict,
    as_content: str,
    outdir: str,
    details_path: str,
    bb_name: str,
    as_name: str,
    chrom_sizes: Optional[str],
    run_tools: bool,
    enzyme: str = "ascas12a",
    validate: bool = True,
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
    # Same invariant check the chunked merge runs, so a single-pass build cannot ship
    # something a chunked build would have caught (or vice versa).
    if validate:
        stats = validate_track_bed(bed_path, autosql_field_names(as_path), enzyme)
        print(
            f"\n\tinvariants OK: {stats['guides']} guides, "
            f"{stats['distinct_percentiles']} distinct raw->percentile pairs\n",
            flush=True,
        )
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
    enzyme: str = "ascas12a",
    validate: bool = True,
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
    any of the enzyme's on-target columns (see :data:`ENZYME_TRACK_SPECS`).
    ``start``/``stop`` bound the 23-nt spacer
    (0-based, + strand). Returns a dict of the artifacts written.

    ``genome_fa`` (an indexed FASTA, ``.fai`` alongside) enables the AsCas12a
    TTTV correction: crispr-ots scores the AsCas12a (2xNLS) CFD PAM-blind, so
    without the reference PAM a TTTT off-target (which AsCas12a cannot cleave)
    would wrongly count against the AsCas12a specificity. With the genome, TTTT
    off-target CFDs are zeroed in the details popup (Layer 2) and the AsCas12a
    close-off-target grey is driven by TTTV-restricted mismatch counts (Layer 3).
    """
    # Per-enzyme filenames + on-target columns. Every track is single-enzyme: one CFD
    # matrix, one off-target index, one specificity axis. Which cleavage models it carries
    # is the only per-enzyme difference (see ENZYME_TRACK_SPECS).
    if enzyme not in ENZYME_TRACK_SPECS:
        raise ValueError(
            f"unknown enzyme {enzyme!r}; expected one of {sorted(ENZYME_TRACK_SPECS)}. "
            "Every track is single-enzyme."
        )
    bb_name = bb_name or f"{enzyme}.bb"
    details_name = details_name or f"{enzyme}_crisprDetails.tab"
    as_name = as_name or f"{enzyme}Targets.as"
    score_cols = ordered_ontarget_cols(enzyme, guide_df.columns)

    os.makedirs(outdir, exist_ok=True)
    mode1 = pd.read_csv(enum_prefix)
    sidecar_path = enum_prefix + ".guides.tsv"
    mode2_path = enum_prefix + ".ot.tsv"

    # --- crisprDetails.tab + offsets ---
    # The ascas12a path uses a TTTV index (no TTTT off-targets exist), so no genome
    # is needed and Mode-1 counts are already TTTV. enCas12a scans the broad index and
    # needs the genome to read each off-target's reference PAM.
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
    offsets = _build_details(
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

    # --- AsCas12a single-enzyme track (TTTV off-target search) ---
    # One matrix (cfd-cas12a:2xnls), one specificity (TTTV == TTTN under a TTTV index).
    # Two AsCas12a on-target models are shown: DeepCpf1 (Kim 2018) and the AsCas12a
    # seq-DeepCpf1 variant (Chen 2025). No enCas12a, no TTTN references, no genome.
    if enzyme == "ascas12a":
        spec = pd.to_numeric(g.get("specificity_tttv"), errors="coerce")
        score = np.clip(np.rint(spec.fillna(0).to_numpy() * 100), 0, 1000).astype(int)
        # Mode-1 mismatch_counts are TTTV because the INDEX is TTTV, not because the
        # scanner restricts them: `mm_counts` is incremented for every site it is handed.
        # `crot_index.assert_enzyme_pams` is what makes that true, and the close-off-target
        # grey reads them directly with no genome re-derivation. Scanned against a TTTN
        # index these counts include TTTT sites AsCas12a cannot cleave, which greys good
        # guides: 8.2% of hg38 chr1 guides, measured 2026-10-22.
        buckets = _mismatch_buckets(g["mismatch_counts"])
        mm0, mm1, mm2 = buckets["mm0"].to_numpy(), buckets["mm1"].to_numpy(), buckets["mm2"].to_numpy()
        ot_cols = ordered_ontarget_cols("ascas12a", g.columns)
        colours, display, hover = _ontarget_block(
            g, "ascas12a", ot_cols, dropped, mm0, mm1, mm2, spec.to_numpy(), score
        )

        out = {
            "#chr": g["chrom"].to_numpy(),
            "start": chrom_start,
            "stop": chrom_end,
            "name": "",
            "score": score,
            "strand": g["strand"].to_numpy(),
            "thickStart": thick_start,
            "thickEnd": thick_end,
            **{f: colours[f] for f in COLOR_FIELDS},
            "guideSeq": g["guideSeq"].to_numpy(),
            "pam": g["pam"].to_numpy(),
            "unique": unique,
            "AsCas12a_specificity": _fmt_raw(spec).to_numpy(),
            **display,
        }
        out["_mouseOver"] = hover
        out["_offset"] = [offsets.get(gid, 0) for gid in g["id"].to_numpy()]
        as_content = build_autosql_ascas12a(ot_cols)
        return _write_track(
            out, as_content, outdir, details_path, bb_name, as_name, chrom_sizes,
            run_tools, enzyme=enzyme, validate=validate,
        )

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

        # Same universal block as AsCas12a: three colourings, split raw/...Perc columns,
        # raw-x100 specificity in the hover, percentiles on the displayed value so the
        # merge can re-rank them genome-wide. The only enzyme-specific inputs are which
        # cleavage models were scored and which of the four specificities drives score.
        ot_cols = ordered_ontarget_cols("encas12a", g.columns)
        colours, display, hover = _ontarget_block(
            g, "encas12a", ot_cols, nonuniq, mm0, mm1, mm2, score_spec, score
        )

        out = {
            "#chr": g["chrom"].to_numpy(),
            "start": chrom_start,
            "stop": chrom_end,
            "name": "",
            "score": score,
            "strand": g["strand"].to_numpy(),
            "thickStart": thick_start,
            "thickEnd": thick_end,
            **{f: colours[f] for f in COLOR_FIELDS},
            "guideSeq": g["guideSeq"].to_numpy(),
            "pam": g["pam"].to_numpy(),
            "unique": unique_en,
            "enCas12a_spec_fa06": _fmt_raw(pd.Series(s06)).to_numpy(),
            "enCas12a_spec_fa06_wt": _fmt_raw(pd.Series(s06w)).to_numpy(),
            "enCas12a_spec_fa03": _fmt_raw(pd.Series(s03)).to_numpy(),
            "enCas12a_spec_fa03_wt": _fmt_raw(pd.Series(s03w)).to_numpy(),
            **display,
        }
        out["_mouseOver"] = hover
        out["_offset"] = [offsets.get(gid, 0) for gid in g["id"].to_numpy()]
        as_content = build_autosql_encas12a(ot_cols)
        return _write_track(
            out, as_content, outdir, details_path, bb_name, as_name, chrom_sizes,
            run_tools, enzyme=enzyme, validate=validate,
        )



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
# Greys a re-grade must PRESERVE. They encode specificity and uniqueness, which
# _item_rgb applies after the ramp and which are genome-wide to begin with, so a cell
# holding one of these is not a ramp colour and must not be recomputed.
#
# `90,90,90` is deliberately NOT here. It was the old marginal (0.20-0.50 specificity)
# tier, which no longer exists, so re-grading an older track converts those cells to
# the ramp colour the guide should now show. That is sound because the later overrides
# win: a cell still reading `90,90,90` means the guide is unique and has no close
# off-target, so the pure ramp is exactly what belongs there.
_RAMP_GREYS = frozenset({b"150,150,150", b"120,120,120", b"50,50,50"})


def track_regrade_spec(enzyme: str, as_fields: Sequence[str]) -> Optional[dict]:
    """Resolve the genome-wide percentile re-grade rule for ``enzyme`` against the
    fields a built track actually declares.

    The rule is derived, not tabulated, so it cannot drift from what
    :func:`build_track` wrote. :data:`ENZYME_TRACK_SPECS` lists each enzyme's
    on-target models in display order; this keeps the ones whose column is present
    in ``as_fields`` and pairs them positionally with :data:`EFF_COLOR_FIELDS`. A run
    that scored one model, or three, therefore re-grades exactly what it has: extra
    cleavage scores are a per-enzyme difference in *which* models exist, never a
    difference in how they are graded or displayed.

    Returns ``None`` for an enzyme with no entry, which :func:`build_track` rejects
    outright; the merge keeps the check as a guard rather than silently not re-grading.
    """
    espec = ENZYME_TRACK_SPECS.get(enzyme)
    if espec is None:
        return None
    models = []
    for col, label, dec in espec["ontarget"]:
        base = _cas9_style_name(col)
        if base not in as_fields:
            continue
        colour = EFF_COLOR_FIELDS[len(models)] if len(models) < len(EFF_COLOR_FIELDS) else None
        if colour is None or colour not in as_fields:
            continue
        models.append((base, label, colour, dec))
    if not models:
        return None
    return {
        "unique": "unique",
        # The hover's specificity is the bed9 ``score`` column read straight back
        # (raw specificity x100), NOT a rank, so it cannot drift from column 5. Only
        # the efficiency models below are ranked and therefore need re-grading.
        "spec_score": "score",
        "spec_label": espec["spec_label"],
        "score_spec_col": espec["score_spec"],
        "hover": "_mouseOver",
        # (raw field, hover label, ramp colour field, displayed decimals)
        "models": models,
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
    O(distinct) in memory, and since each model is written at its declared precision
    (0 or 2 decimals) a genome of 10^8 guides has only ~10^2 distinct keys.

    This reads the WRITTEN bytes, which is why re-ranking on the rounded display
    needs no change here: the chunks already wrote the rounded value, so the
    genome-wide rank the merge computes is the rank of what the browser shows."""
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


# Every colour a track is allowed to emit. Anything else in a colour column means
# the colour logic and this list have drifted apart, which is what shipped the
# retired "90,90,90" marginal-specificity grey into five genomes before anyone looked.
def _is_fixed_decimals(v: bytes, decimals: int) -> bool:
    """Is ``v`` exactly what :func:`_fmt_score` writes at ``decimals`` places?

    Pure byte inspection rather than a ``float()`` plus reformat, because the
    validator touches every guide in the genome and a per-row round trip would cost
    more than the rest of the pass put together."""
    if v[:1] == b"-":
        v = v[1:]
    if decimals == 0:
        return v.isdigit()
    i = v.find(b".")
    return i > 0 and len(v) - i - 1 == decimals and v[:i].isdigit() and v[i + 1 :].isdigit()


RAMP_COLOURS = frozenset({b"255,100,100", b"255,255,0", b"0,200,0"})
VALID_COLOURS = RAMP_COLOURS | _RAMP_GREYS
VALID_SPEC_COLOURS = frozenset({b"255,0,0", b"128,128,0", b"0,255,0"})


def validate_track_bed(bed_path: str, as_fields: Sequence[str], enzyme: str) -> Dict[str, int]:
    """Check a finished track BED against the invariants the track is meant to hold,
    and raise :class:`ValueError` listing every violation found.

    This runs inside the build, on the final sorted BED, so a track that violates an
    invariant fails the run instead of being discovered later and patched by hand. It
    is a single streaming pass and keeps only per-value tallies, so it costs one read
    of the BED regardless of genome size.

    What it asserts:

    * every row has exactly the field count the autoSql declares;
    * the three colour columns hold only colours the scheme can produce, so a retired
      or stray colour cannot ship;
    * a colour column is a ramp colour only where the guide is unique and free of
      close off-targets, since the greys override the ramp;
    * each raw cleavage score maps to exactly ONE percentile, which is what makes the
      percentile genome-wide rather than chunk-local;
    * each cleavage score is written to exactly the decimals its model declares in
      :data:`ENZYME_TRACK_SPECS`. Paired with the check above this pins the whole
      display contract: the percentile must be the rank of the number actually shown,
      so a column written at one precision and ranked at another fails here;
    * a percentile is blank exactly when its raw value is blank;
    * ``score`` agrees with the specificity column to within 1, so column 5 and the
      specificity cannot drift apart. The tolerance is not slack: ``score`` is rounded
      from the full-precision float while the column is written to 4 decimals, so when
      ``displayed x 100`` lands exactly on a half the two legitimately round opposite
      ways and differ by 1. Measured on sacCer3, 113 of 590,425 guides sit on such a
      half (e.g. a displayed 0.905 whose score is 91). Anything beyond 1 is real drift;
    * the hover reports that same ``score``, not a percentile.
    """
    spec = track_regrade_spec(enzyme, as_fields)
    if spec is None:
        raise ValueError(f"cannot validate enzyme {enzyme!r}: no track spec")
    idx = {f: i for i, f in enumerate(as_fields)}
    n_f = len(as_fields)
    i_spec_col = idx["_specColor"]
    i_uniq, i_hover, i_score = idx[spec["unique"]], idx[spec["hover"]], idx["score"]
    spec_raw_col = spec["score_spec_col"]
    if spec_raw_col not in idx:
        raise ValueError(
            f"{enzyme} track is missing its score-driving specificity column {spec_raw_col!r}"
        )
    i_spec_raw = idx[spec_raw_col]
    label = spec["spec_label"].encode() + b": "
    models = [
        (idx[raw], idx.get(raw + "Perc"), idx[col], raw, dec)
        for raw, _lab, col, dec in spec["models"]
    ]

    bad_colour: Dict[bytes, int] = {}
    pct_of: Dict[str, Dict[bytes, bytes]] = {raw: {} for _ir, _ip, _ic, raw, _d in models}
    multi: Dict[str, set] = {raw: set() for _ir, _ip, _ic, raw, _d in models}
    bad_fmt: Dict[str, int] = {raw: 0 for _ir, _ip, _ic, raw, _d in models}
    n = n_blank_mismatch = n_score_mismatch = n_hover_mismatch = n_ramp_on_grey = 0

    with open(bed_path, "rb") as fh:
        for row in fh:
            row = row.rstrip(b"\n")
            if not row:
                continue
            f = row.split(b"\t")
            if len(f) != n_f:
                raise ValueError(
                    f"{bed_path}: row {n} has {len(f)} fields, autoSql declares {n_f}"
                )
            n += 1
            greyed = f[idx[spec["unique"]]] == b"False"
            for c in (f[i] for _ir, _ip, i, _raw, _d in models):
                if c not in VALID_COLOURS:
                    bad_colour[c] = bad_colour.get(c, 0) + 1
                elif greyed and c in RAMP_COLOURS:
                    n_ramp_on_grey += 1
            if f[i_spec_col] not in VALID_SPEC_COLOURS:
                bad_colour[f[i_spec_col]] = bad_colour.get(f[i_spec_col], 0) + 1
            for i_raw, i_pct, _i_col, raw, dec in models:
                if f[i_raw] and not _is_fixed_decimals(f[i_raw], dec):
                    bad_fmt[raw] += 1
                if i_pct is None:
                    continue
                rv, pv = f[i_raw], f[i_pct]
                if (rv == b"") != (pv == b""):
                    n_blank_mismatch += 1
                    continue
                if rv == b"":
                    continue
                seen = pct_of[raw].setdefault(rv, pv)
                if seen != pv:
                    multi[raw].add(rv)
            # score vs the specificity column, to within the 4-decimal rounding (see
            # the docstring), and the hover repeats score verbatim
            if not greyed and f[i_spec_raw]:
                if abs(int(f[i_score]) - float(f[i_spec_raw]) * 100.0) > 1.0:
                    n_score_mismatch += 1
            if not greyed and not f[i_hover].startswith(label + f[i_score]):
                n_hover_mismatch += 1

    problems = []
    if bad_colour:
        shown = ", ".join(f"{c.decode()} x{k}" for c, k in sorted(bad_colour.items())[:6])
        problems.append(f"{sum(bad_colour.values())} cells hold a colour the scheme cannot produce ({shown})")
    if n_ramp_on_grey:
        problems.append(f"{n_ramp_on_grey} non-unique guides carry a ramp colour instead of a grey")
    for raw, vals in multi.items():
        if vals:
            ex = ", ".join(v.decode() for v in sorted(vals)[:3])
            problems.append(
                f"{raw}: {len(vals)} raw values map to more than one percentile "
                f"(e.g. {ex}); percentiles are not genome-wide"
            )
    declared = {raw: dec for _ir, _ip, _ic, raw, dec in models}
    for raw, k in bad_fmt.items():
        if k:
            problems.append(
                f"{raw}: {k} values are not written to the {declared[raw]} decimals "
                f"ENZYME_TRACK_SPECS declares for that model"
            )
    if n_blank_mismatch:
        problems.append(f"{n_blank_mismatch} rows have a raw value without a percentile, or the reverse")
    if n_score_mismatch:
        problems.append(
            f"{n_score_mismatch} rows where score and the specificity column differ by more than 1"
        )
    if n_hover_mismatch:
        problems.append(f"{n_hover_mismatch} rows where the hover specificity does not equal score")
    if problems:
        raise ValueError(
            f"track invariants violated in {bed_path} ({n} guides):\n\t"
            + "\n\t".join(problems)
        )
    return {
        "guides": n,
        "distinct_percentiles": sum(len(v) for v in pct_of.values()),
    }


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
    i_spec = idx[spec["spec_score"]]
    i_uniq = idx[spec["unique"]]
    i_hover = idx[spec["hover"]]
    label_spec = spec["spec_label"].encode()
    models = [
        (idx[raw], idx.get(raw + "Perc"), idx[col], label.encode(), lut.get(raw, {}))
        for raw, label, col, _dec in spec["models"]
    ]

    n = 0
    with open(bed_in, "rb") as fin, open(bed_out, "wb") as fout:
        for row in fin:
            row = row.rstrip(b"\n")
            if not row:
                continue
            f = row.split(b"\t")
            # hover: "Not unique" for a non-unique guide, else the raw specificity
            # x100 read straight from the score column, then one entry per model that
            # scored. Specificity is absolute, not ranked, so nothing to re-grade there.
            if f[i_uniq] == b"False":
                parts = [b"Not unique"]
            else:
                parts = [label_spec + b": " + f[i_spec]]
            for i_raw, i_pct, i_col, label, model_lut in models:
                pct = model_lut.get(f[i_raw]) if f[i_raw] else None
                if i_pct is not None:
                    f[i_pct] = b"" if pct is None else b"%d" % round(pct)
                if f[i_col] not in _RAMP_GREYS:
                    f[i_col] = _ramp_colour(pct)
                if f[i_raw]:
                    # "%" because the figure IS a percentile, as on the Cas9 track
                    parts.append(label + b": " + b"%d%%" % round(pct))
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
    details_name: str = "ascas12a_crisprDetails.tab",
    bb_name: str = "ascas12a.bb",
    as_name: str = "ascas12aTargets.as",
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
    entry in :data:`ENZYME_TRACK_SPECS` is merged without re-grading and warns, because
    its percentiles would stay chunk-local.

    Returns the artifact paths.
    """
    header = b"_mismatchCounts\t_crisprOfftargets\n"

    # Columns whose percentiles the merge must re-rank genome-wide, and the raw
    # columns to histogram while streaming the per-chunk BEDs.
    as_fields = autosql_field_names(as_file)
    spec = track_regrade_spec(enzyme, as_fields)
    if spec is not None:
        needed = [spec["spec_score"], spec["unique"], spec["hover"]]
        needed += [raw for raw, _lab, _col, _dec in spec["models"]]
        needed += [col for _raw, _lab, col, _dec in spec["models"]]
        missing = [f for f in needed if f not in as_fields]
        if missing:
            raise ValueError(f"{as_file} is missing fields needed to re-grade percentiles: {missing}")
    else:
        print(
            f"\n\tWARNING: no percentile re-grading rule for enzyme={enzyme!r}; the merged track's "
            f"percentiles and efficiency-ramp colours stay ranked WITHIN each chunk, not genome-wide.\n",
            flush=True,
        )
    # only the ranked columns need a histogram; specificity is absolute
    # The chunk BEDs already store each raw column AT its displayed precision, so
    # histogramming the written bytes ranks the displayed value, which is what the
    # unchunked path ranks too. That is why no decimals are needed here.
    count_fields = [raw for raw, _lab, _col, _dec in spec["models"]] if spec else []
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
    stats = validate_track_bed(sorted_bed, as_fields, enzyme)
    print(
        f"  invariants OK: {stats['guides']} guides, "
        f"{stats['distinct_percentiles']} distinct raw->percentile pairs",
        flush=True,
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
