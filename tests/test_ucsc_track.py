"""Unit tests for crisprware.ucsc_track (UCSC Cas12a track assembly).

Exercises the adapter on synthetic crispr-ots streaming output: the strand-aware
5'-PAM bigBed geometry, the off-target browser-pos mapping (+: pos=start,
-: pos=start+3), _mismatchCounts/_crisprOfftargets assembly, and the _offset ->
byte-offset round-trip. No TensorFlow / GPU / UCSC tools required (run_tools=False).
"""

import os
import re
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from crisprware import ucsc_track  # noqa: E402



def _as_index(as_path):
    """Field name -> BED column index, read from the track's own autoSql.

    Indexing by name rather than position keeps these tests from breaking every
    time a column is added, which is what happened when the two alternative
    colour fields became universal."""
    return {f: i for i, f in enumerate(ucsc_track.autosql_field_names(as_path))}

def _write_enum(prefix):
    """Two guides: g_uniq (unique, off-targets at + and -), g_dup (perfect-match
    duplicate -> dropped)."""
    # Mode-1 aggregated
    with open(prefix, "w") as f:
        f.write("id,specificity_tttn,specificity_tttv,max_cfd,off_target_count,mismatch_counts,saturated,dropped\n")
        f.write("g_uniq,0.5,0.6,0.9,3,0;1;0;2;0,0,0\n")
        f.write("g_dup,,,,,,0,1\n")
    # sidecar (non-dropped only): guide_id 0 -> g_uniq
    with open(prefix + ".guides.tsv", "w") as f:
        f.write("guide_id\tid\tchrom\tposition\tstrand\tsequence\ton_count\n")
        f.write("0\tg_uniq\tchr22\t1000\t+\tACGTACGTACGTACGTACGTACG\t1\n")
    # Mode-2 per-off-target: one + (start 500) and one - (start 800)
    with open(prefix + ".ot.tsv", "w") as f:
        f.write("guide_id\tchrom\tstart\tstrand\tmm\tcfd\n")
        f.write("0\tchr1\t500\t+\t1\t0.80\n")
        f.write("0\tchr3\t800\t-\t3\t0.30\n")


def _guide_df():
    # start/stop bound the 23-nt spacer on the + strand (as generate_guides emits)
    return pd.DataFrame(
        {
            "id": ["g_uniq", "g_dup"],
            "chrom": ["chr22", "chr22"],
            "start": [104, 204],  # spacer start
            "stop": [127, 227],  # spacer start + 23
            "strand": ["+", "-"],
            "guideSeq": ["ACGTACGTACGTACGTACGTACG", "TTTTGGGGCCCCAAAATTTTGGG"],
            "pam": ["TTTA", "TTTC"],
            "deepcpf1_score": [42.0, 10.0],
            "enseq_deepcpf1_score": [0.6, 0.2],
        }
    )


def test_track_geometry_and_offsets(tmp_path):
    prefix = str(tmp_path / "off.csv")
    _write_enum(prefix)
    art = ucsc_track.build_track(_guide_df(), prefix, str(tmp_path), run_tools=False, enzyme="ascas12a")

    bed = [line.rstrip("\n").split("\t") for line in open(art["bed"])]
    assert len(bed) == 2
    f = _as_index(art["autosql"])
    by_chromstart = {int(r[1]): r for r in bed}

    # g_uniq (+): spacer [104,127) -> site [100,127) (27nt), thick [104,127) (23nt)
    plus = by_chromstart[100]
    assert int(plus[2]) - int(plus[1]) == 27  # chromEnd - chromStart
    assert int(plus[7]) - int(plus[6]) == 23  # thickEnd - thickStart
    assert plus[6] == "104" and plus[7] == "127"  # thick = spacer
    assert plus[f["pam"]] == "TTTA"  # 4-nt PAM
    assert plus[5] == "+"
    assert int(plus[-1]) > 0  # _offset set (has off-targets)

    # g_dup (-): spacer [204,227) -> site [204,231) (27nt), thick [204,227) (PAM on the right)
    minus = next(r for r in bed if r[5] == "-")
    assert int(minus[1]) == 204 and int(minus[2]) == 231  # chromStart, chromEnd
    assert int(minus[6]) == 204 and int(minus[7]) == 227  # thick = spacer
    assert minus[f["itemRgb"]] == "150,150,150"  # grey (non-unique/dropped)
    assert minus[-1] == "0"  # dropped -> no crisprDetails line

    # All geometry invariants
    for r in bed:
        assert int(r[2]) - int(r[1]) == 27
        assert int(r[7]) - int(r[6]) == 23
        assert len(r[f["pam"]]) == 4


def test_offtarget_browser_pos_and_offset_seek(tmp_path):
    prefix = str(tmp_path / "off.csv")
    _write_enum(prefix)
    art = ucsc_track.build_track(_guide_df(), prefix, str(tmp_path), run_tools=False, enzyme="ascas12a")

    # find g_uniq's _offset
    bed = [line.rstrip("\n").split("\t") for line in open(art["bed"])]
    off = next(int(r[-1]) for r in bed if r[1] == "100")

    raw = open(art["details_tab"], "rb").read()
    # header at byte 0; the recorded offset points at g_uniq's line
    line = raw[off:].split(b"\n", 1)[0].decode()
    counts, lst = line.split("\t")
    # The histogram counts the off-targets that are actually listed, NOT Mode-1's
    # totals. Mode-1 reports 0;1;0;2;0 (3 hits) but one of the two mm=3 hits is
    # below the scan's CFD floor and never reaches Mode-2, so it is neither shown
    # nor counted (Cas9/CRISPOR convention).
    assert counts == "0,1,0,1,0"
    # off-targets sorted by score desc; pos = crispr-ots start for BOTH strands
    # (the full 27-nt site's leftmost base): + at 500, - at 800.
    entries = lst.split("|")
    assert entries[0] == "chr1;500+;800"  # cfd 0.80 -> 800; + strand pos=start
    assert entries[1] == "chr3;800-;300"  # cfd 0.30 -> 300; - strand pos=start (was start+3)
    # The invariant the browser relies on: every counted off-target is shown.
    assert sum(int(x) for x in counts.split(",")) == len(entries)


def test_blank_list_when_too_many(tmp_path):
    prefix = str(tmp_path / "off.csv")
    _write_enum(prefix)
    # off_target_count 3 > blank_threshold 1 -> list blanked. The histogram keeps
    # the true above-floor totals (2), which is what makes the browser print its
    # "too many off-targets" message instead of showing a truncated list.
    art = ucsc_track.build_track(_guide_df(), prefix, str(tmp_path), blank_threshold=1, run_tools=False, enzyme="ascas12a")
    raw = open(art["details_tab"], "rb").read().decode()
    body = raw.split("\n", 1)[1]
    counts, lst = body.split("\n")[0].split("\t")
    assert counts == "0,1,0,1,0"
    assert lst == ""  # blanked


def test_over_cap_emits_empty_list_with_full_histogram(tmp_path):
    """Over the per-guide cap, follow Cas9: empty list, true above-floor counts.

    A truncated list would silently contradict the histogram; an empty one tells
    the browser to show the "use the Crispor.org link" message.
    """
    prefix = str(tmp_path / "off.csv")
    _write_enum(prefix)
    art = ucsc_track.build_track(_guide_df(), prefix, str(tmp_path), list_cap=1, run_tools=False, enzyme="ascas12a")
    body = open(art["details_tab"]).read().split("\n", 1)[1]
    counts, lst = body.split("\n")[0].split("\t")
    assert lst == ""  # 2 off-targets > cap of 1 -> not truncated to 1, emptied
    assert counts == "0,1,0,1,0"  # still the true above-floor histogram

    # At a cap that fits, the list is shown and the invariant holds exactly.
    art2 = ucsc_track.build_track(_guide_df(), prefix, str(tmp_path / "b"), list_cap=2, run_tools=False, enzyme="ascas12a")
    counts2, lst2 = open(art2["details_tab"]).read().split("\n", 1)[1].split("\n")[0].split("\t")
    assert sum(int(x) for x in counts2.split(",")) == len(lst2.split("|")) == 2


def test_all_offtargets_below_floor_report_zero(tmp_path):
    """A guide whose hits are all sub-floor counts zero, not Mode-1's total."""
    prefix = str(tmp_path / "off.csv")
    _write_enum(prefix)
    open(prefix + ".ot.tsv", "w").write("guide_id\tchrom\tstart\tstrand\tmm\tcfd\n")  # no Mode-2 rows
    art = ucsc_track.build_track(_guide_df(), prefix, str(tmp_path), run_tools=False, enzyme="ascas12a")
    body = open(art["details_tab"]).read().split("\n", 1)[1]
    counts, lst = body.split("\n")[0].split("\t")
    assert lst == ""
    assert counts == "0,0,0,0,0"  # not Mode-1's "0,1,0,2,0"


def _write_fasta(path, contigs):
    """Single-line-per-contig FASTA + a matching .fai (offsets read from the bytes
    actually written), so Genome works without samtools in the test env."""
    content = "".join(f">{name}\n{seq}\n" for name, seq in contigs)
    with open(path, "w") as f:
        f.write(content)
    b = content.encode()
    with open(path + ".fai", "w") as f:
        for name, seq in contigs:
            seq_off = b.index(f">{name}\n".encode()) + len(name) + 2  # '>' + name + '\n'
            f.write(f"{name}\t{len(seq)}\t{seq_off}\t{len(seq)}\t{len(seq) + 1}\n")


def _write_ascas12a_enum(prefix):
    """AsCas12a single-matrix scan against a TTTV-only index: every off-target is
    TTTV (single ``cfd`` column, no ``cfd2``), Mode-1 has no ``_2xnls`` columns.
    g_hi: unique, one distant (3-mm) off-target. g_close: unique, one 1-mm close
    off-target. g_dup: dropped (non-unique)."""
    with open(prefix, "w") as f:
        f.write("id,specificity_tttn,specificity_tttv,max_cfd,off_target_count,mismatch_counts,saturated,dropped\n")
        f.write("g_hi,0.9,0.9,0.4,1,0;0;0;1;0,0,0\n")
        f.write("g_close,0.6,0.6,0.5,1,0;1;0;0;0,0,0\n")
        f.write("g_dup,,,,,,0,1\n")
    with open(prefix + ".guides.tsv", "w") as f:
        f.write("guide_id\tid\tchrom\tposition\tstrand\tsequence\ton_count\n")
        f.write("0\tg_hi\tchr22\t1000\t+\tACGTACGTACGTACGTACGTACG\t1\n")
        f.write("1\tg_close\tchr22\t2000\t+\tGGGGCCCCAAAATTTTACGTACG\t1\n")
    with open(prefix + ".ot.tsv", "w") as f:
        f.write("guide_id\tchrom\tstart\tstrand\tmm\tcfd\n")
        f.write("0\tchr1\t500\t+\t3\t0.40\n")
        f.write("1\tchr3\t800\t-\t1\t0.50\n")


def _ascas12a_guide_df():
    return pd.DataFrame(
        {
            "id": ["g_hi", "g_close", "g_dup"],
            "chrom": ["chr22", "chr22", "chr22"],
            "start": [104, 204, 304],
            "stop": [127, 227, 327],
            "strand": ["+", "+", "-"],
            "guideSeq": [
                "ACGTACGTACGTACGTACGTACG",
                "GGGGCCCCAAAATTTTACGTACG",
                "TTTTGGGGCCCCAAAATTTTGGG",
            ],
            "pam": ["TTTA", "TTTC", "TTTG"],
            "deepcpf1_score": [30.0, 55.0, 5.0],
            "ascas12a_deepcpf1_score": [80.0, 45.0, 10.0],
        }
    )


def test_build_track_ascas12a_single_enzyme(tmp_path):
    prefix = str(tmp_path / "off.csv")
    _write_ascas12a_enum(prefix)
    art = ucsc_track.build_track(
        _ascas12a_guide_df(), prefix, str(tmp_path), run_tools=False, enzyme="ascas12a"
    )

    # per-enzyme filenames
    assert art["autosql"].endswith("ascas12aTargets.as")
    assert art["details_tab"].endswith("ascas12a_crisprDetails.tab")

    # schema: single AsCas12a specificity, NO TTTN / enCas12a anywhere
    aso = open(art["autosql"]).read()
    assert "AsCas12a_specificity" in aso
    # BOTH AsCas12a on-target models, each as "label|description" so hgc renders a
    # short field label with the citation as a note beneath it
    assert "Efficiency: Kim et al. 2018 Score (raw)" in aso
    assert "Efficiency: Kim et al. 2018 Score (percentile)" in aso
    assert "Efficiency: Chen et al. 2025 Score (raw)" in aso
    assert "Efficiency: Chen et al. 2025 Score (percentile)" in aso
    assert "AsCas12a aggregate CFD specificity" in aso
    # the obvious name is wrong: what runs is Seq-deepCpf1, so nothing claims "DeepCpf1"
    assert "DeepCpf1 (Kim" not in aso
    for banned in ("TTTN", "EnCas12a", "enCas12a", "unique_TTTV", "unique_TTTN", "specificity_tttv_"):
        assert banned not in aso

    bed = [line.rstrip("\n").split("\t") for line in open(art["bed"])]
    # bed9 + _crisprDeepCpf1Color,_specColor,guideSeq,pam,unique,spec,
    # deepcpf1,deepcpf1Perc,ascas12aDeepcpf1,ascas12aDeepcpf1Perc,_mouseOver,_offset
    assert all(len(r) == 21 for r in bed)
    assert all(r[3] == "" for r in bed)  # name is empty (Cas9-matching), not a click-fragile " "
    assert all("%" not in r[14] for r in bed)  # AsCas12a_specificity is raw CFD (0-1), not a degenerate percentile
    # raw and percentile are SEPARATE numeric fields (Cas9 fusi/fusiPerc), not one string
    for r in bed:
        for idx in (15, 16, 17, 18):
            assert "%" not in r[idx] and "(" not in r[idx], r[idx]
    by_score = {int(r[4]): r for r in bed}
    assert set(by_score) == {90, 60, 0}  # score = round(specificity_tttv * 100)

    # itemRgb (r[8]) keeps the Cas9 default blended scheme; _specColor (r[9]) is the
    # additional specificity-only coloring.
    hi, close, dup = by_score[90], by_score[60], by_score[0]
    assert hi[8] in {"255,100,100", "255,255,0", "0,200,0"}  # spec 0.9, no close OT -> efficiency color
    assert close[8] == "120,120,120"  # 1-mm close off-target (Mode-1 TTTV counts) -> grey
    assert dup[8] == "150,150,150"  # dropped -> grey
    # _crisprDeepCpf1Color (r[9]) uses the same scheme on the second on-target model
    assert close[9] == "120,120,120" and dup[9] == "150,150,150"
    # _specColor (r[10]) is specificity ALONE: Cas9 keeps a colour even when greyed
    assert hi[10] == "0,255,0"  # spec 0.90 -> >70 -> green
    assert dup[10] == "255,0,0"  # non-unique -> worst tier, still coloured
    assert dup[-1] == "0"  # dropped -> no details line

    hovers = [r[-2] for r in bed]  # _mouseOver
    assert any("AsCas12a Spec" in h for h in hovers)
    assert all("EnCas12a" not in h and "TTTN" not in h for h in hovers)
    assert "Not unique" in dup[-2]


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_encas12a_enum(prefix, cfd06=0.10, cfd03=0.05):
    """enCas12a single-matrix scan against the broad fa>0.3 index: one cfd column, no cfd2.
    g_uniq: unique, two DISTANT (3-mm) off-targets -- chr1:500(+) at a high-activity (fa>0.6)
    PAM, chr3:800(-) at a moderate (fa>0.3-only) PAM. g_dup: dropped (non-unique)."""
    with open(prefix, "w") as f:
        f.write("id,specificity_tttn,specificity_tttv,max_cfd,off_target_count,mismatch_counts,saturated,dropped\n")
        f.write(f"g_uniq,0.9,0.9,{max(cfd06, cfd03)},2,0;0;0;2;0,0,0\n")  # two 3-mm -> no close grey
        f.write("g_dup,,,,,,0,1\n")
    with open(prefix + ".guides.tsv", "w") as f:
        f.write("guide_id\tid\tchrom\tposition\tstrand\tsequence\ton_count\n")
        f.write("0\tg_uniq\tchr22\t1000\t+\tACGTACGTACGTACGTACGTACG\t1\n")
    with open(prefix + ".ot.tsv", "w") as f:
        f.write("guide_id\tchrom\tstart\tstrand\tmm\tcfd\n")
        f.write(f"0\tchr1\t500\t+\t3\t{cfd06:.2f}\n")
        f.write(f"0\tchr3\t800\t-\t3\t{cfd03:.2f}\n")


def test_build_track_encas12a_four_specificities(tmp_path):
    """enCas12a single-enzyme track: four guidescan specificities (raw + PAM-weighted, at the
    fa>0.6 and fa>0.3 cutoffs) folded in from the Mode-2 stream; col5 = spec_fa0.3_w x100;
    color driven by EnPAM-GB; schema carries no AsCas12a / TTTV / TTTN references."""
    weights = ucsc_track.load_encas12a_weights()
    pam06 = next(p for p, (w, i6, i3) in weights.items() if i6 and i3)  # high-activity, fa>0.6
    pam03 = next(p for p, (w, i6, i3) in weights.items() if i3 and not i6)  # moderate, fa>0.3 only
    w06, w03 = weights[pam06][0], weights[pam03][0]
    cfd06, cfd03 = 0.10, 0.05

    prefix = str(tmp_path / "off.csv")
    _write_encas12a_enum(prefix, cfd06, cfd03)
    # genome: chr1[500:504] = pam06 (+ off-target PAM); chr3[823:827] = revcomp(pam03)
    # (- off-target PAM = revcomp of [pos+23, pos+27)).
    chr1 = "A" * 500 + pam06 + "A" * 40
    chr3 = "A" * 823 + _rc(pam03) + "A" * 13
    fa = str(tmp_path / "g.fa")
    _write_fasta(fa, [("chr1", chr1), ("chr3", chr3)])

    gdf = pd.DataFrame(
        {
            "id": ["g_uniq", "g_dup"],
            "chrom": ["chr22", "chr22"],
            "start": [104, 204],
            "stop": [127, 227],
            "strand": ["+", "-"],
            "guideSeq": ["ACGTACGTACGTACGTACGTACG", "TTTTGGGGCCCCAAAATTTTGGG"],
            "pam": ["TTTA", "TTTC"],
            "enpam_gb_score": [0.90, 0.10],  # g_uniq higher -> pct 100 -> green ramp
            "encas12a_hf1_score": [0.10, 0.80],  # low for g_uniq: proves color != HF1
        }
    )
    art = ucsc_track.build_track(
        gdf, prefix, str(tmp_path), run_tools=False, enzyme="encas12a", genome_fa=fa, encas12a_list_floor=0.0
    )

    # per-enzyme filenames + schema
    assert art["autosql"].endswith("encas12aTargets.as")
    assert art["details_tab"].endswith("encas12a_crisprDetails.tab")
    aso = open(art["autosql"]).read()
    for col in ("enCas12a_spec_fa06", "enCas12a_spec_fa06_wt", "enCas12a_spec_fa03", "enCas12a_spec_fa03_wt"):
        assert col in aso
    assert "EnPAM-GB" in aso and "HF1" in aso
    # no references to the OTHER enzyme's specificity ("enAsCas12a", the enzyme's own name, is fine)
    assert "AsCas12a" not in aso.replace("enAsCas12a", "")
    for banned in ("TTTV", "TTTN", "specificity_tttv"):
        assert banned not in aso

    bed = [line.rstrip("\n").split("\t") for line in open(art["bed"])]
    # bed9 + 2 alt colours + guideSeq,pam,unique + 4 spec + 2 ot x (raw,Perc) + _mouseOver,_offset
    assert all(len(r) == 24 for r in bed)
    assert all(r[3] == "" for r in bed)  # empty name (not a click-fragile " ")
    uniq = next(r for r in bed if r[5] == "+")
    dup = next(r for r in bed if r[5] == "-")

    # four specificities (1 / (on_count + weighted CFD sums)); fa>0.6 counts only the fa>0.6
    # off-target, fa>0.3 counts both.
    oc = 1
    f = _as_index(art["autosql"])
    exp = {
        "enCas12a_spec_fa06": 1 / (oc + cfd06),
        "enCas12a_spec_fa06_wt": 1 / (oc + cfd06 * w06),
        "enCas12a_spec_fa03": 1 / (oc + cfd06 + cfd03),
        "enCas12a_spec_fa03_wt": 1 / (oc + cfd06 * w06 + cfd03 * w03),
    }
    for col, e in exp.items():
        got = uniq[f[col]]
        assert float(got) == pytest.approx(e, abs=1e-3), f"{col}: {got} != {e:.4f}"
    # col5 = spec_fa0.3_w x100
    assert int(uniq[4]) == round(exp["enCas12a_spec_fa03_wt"] * 100)

    # color: no close off-target (3-mm) + spec well above the grey threshold -> EnPAM-GB
    # efficiency ramp; g_uniq has the top EnPAM-GB percentile -> green.
    assert uniq[f["itemRgb"]] == "0,200,0"
    assert int(uniq[-1]) > 0  # _offset set (has off-targets)
    assert "enCas12a Spec. Score (fa>0.3 wt)" in uniq[-2] and "EnPAM-GB" in uniq[-2]

    # raw and percentile are SEPARATE columns, as on the AsCas12a track: no packed
    # "NN% (raw)" string in any DATA column. The hover is excluded because it carries
    # the percentile with a % suffix on purpose, matching the Cas9 track's
    # "Doench 2016: 18%"; forbidding % outright would ban that too.
    assert all("%" not in c for c in uniq[:-2]), "a data column holds a packed percent string"
    assert uniq[f["enpamGbPerc"]].isdigit()
    assert "EnPAM-GB: " + uniq[f["enpamGbPerc"]] + "%" in uniq[-2]
    assert float(uniq[f["enpamGb"]]) == pytest.approx(0.9, abs=1e-4)

    # dropped guide: greyed, blank specificity columns, no details line
    assert dup[f["itemRgb"]] == "150,150,150" and dup[-1] == "0"
    assert all(dup[f[c]] == "" for c in exp)
    assert "Not unique" in dup[-2]


def test_build_track_encas12a_blanks_perfect_repeat(tmp_path):
    """A guide with a perfect (0-mm) off-target elsewhere (mm0>0) but NOT pre-dropped must be
    treated as non-unique: blank all four specificity columns, unique=False, greyed. Otherwise
    a repeat whose off-target sits at a non-fa0.6 PAM would read a misleading spec_fa0.6=1.0."""
    weights = ucsc_track.load_encas12a_weights()
    pam06 = next(p for p, (w, i6, i3) in weights.items() if i6 and i3)
    prefix = str(tmp_path / "off.csv")
    with open(prefix, "w") as f:
        f.write("id,specificity_tttn,specificity_tttv,max_cfd,off_target_count,mismatch_counts,saturated,dropped\n")
        f.write("g_rep,0.5,0.5,1.0,1,1;0;0;0;0,0,0\n")  # one 0-mm off-target, dropped=0
    with open(prefix + ".guides.tsv", "w") as f:
        f.write("guide_id\tid\tchrom\tposition\tstrand\tsequence\ton_count\n")
        f.write("0\tg_rep\tchr22\t1000\t+\tACGTACGTACGTACGTACGTACG\t1\n")
    with open(prefix + ".ot.tsv", "w") as f:
        f.write("guide_id\tchrom\tstart\tstrand\tmm\tcfd\n")
        f.write("0\tchr1\t500\t+\t0\t1.00\n")  # perfect match elsewhere
    chr1 = "A" * 500 + pam06 + "A" * 40
    fa = str(tmp_path / "g.fa")
    _write_fasta(fa, [("chr1", chr1)])
    gdf = pd.DataFrame(
        {"id": ["g_rep"], "chrom": ["chr22"], "start": [104], "stop": [127], "strand": ["+"],
         "guideSeq": ["ACGTACGTACGTACGTACGTACG"], "pam": ["TTTA"],
         "enpam_gb_score": [0.7], "encas12a_hf1_score": [0.3]}
    )
    art = ucsc_track.build_track(
        gdf, prefix, str(tmp_path), run_tools=False, enzyme="encas12a", genome_fa=fa, encas12a_list_floor=0.0
    )
    r = open(art["bed"]).readline().rstrip("\n").split("\t")
    f = _as_index(art["autosql"])
    assert r[f["unique"]] == "False"
    assert all(r[f[c]] == "" for c in ("enCas12a_spec_fa06", "enCas12a_spec_fa06_wt",
                                       "enCas12a_spec_fa03", "enCas12a_spec_fa03_wt"))
    assert r[4] == "0"  # col5 score 0
    assert r[f["itemRgb"]] == "150,150,150"  # non-unique grey
    assert "Not unique" in r[-2]


def test_encas12a_spec_exact_vs_pilot_formula(tmp_path):
    """The folded-in specificity must equal the standalone pilot formula
    1/(on_count + sum CFD[*weight]) to the emitted precision, for a unique guide with a
    unique on-target and several off-targets at mixed PAMs (regression guard for the
    single-command track vs the validated pilot aggregate_specificity)."""
    weights = ucsc_track.load_encas12a_weights()
    p06 = next(p for p, (w, i6, i3) in weights.items() if i6 and i3)
    p03 = next(p for p, (w, i6, i3) in weights.items() if i3 and not i6)
    ots = [("chr1", 500, "+", 0.40, p06), ("chr1", 600, "+", 0.12, p03), ("chr3", 800, "-", 0.22, p06)]

    prefix = str(tmp_path / "off.csv")
    with open(prefix, "w") as f:
        f.write("id,specificity_tttn,specificity_tttv,max_cfd,off_target_count,mismatch_counts,saturated,dropped\n")
        f.write("g,0.9,0.9,0.4,3,0;0;3;0;0,0,0\n")
    with open(prefix + ".guides.tsv", "w") as f:
        f.write("guide_id\tid\tchrom\tposition\tstrand\tsequence\ton_count\n")
        f.write("0\tg\tchr22\t1000\t+\tACGTACGTACGTACGTACGTACG\t1\n")
    with open(prefix + ".ot.tsv", "w") as f:
        f.write("guide_id\tchrom\tstart\tstrand\tmm\tcfd\n")
        for chrom, pos, strand, cfd, _pam in ots:
            f.write(f"0\t{chrom}\t{pos}\t{strand}\t2\t{cfd:.2f}\n")
    # genome placing each off-target PAM
    seqs = {"chr1": ["A"] * 700, "chr3": ["A"] * 900}
    for chrom, pos, strand, _cfd, pam in ots:
        s = seqs[chrom]
        if strand == "+":
            s[pos : pos + 4] = list(pam)
        else:
            s[pos + 23 : pos + 27] = list(_rc(pam))
    fa = str(tmp_path / "g.fa")
    _write_fasta(fa, [(c, "".join(s)) for c, s in seqs.items()])

    gdf = pd.DataFrame(
        {"id": ["g"], "chrom": ["chr22"], "start": [104], "stop": [127], "strand": ["+"],
         "guideSeq": ["ACGTACGTACGTACGTACGTACG"], "pam": ["TTTA"],
         "enpam_gb_score": [0.5], "encas12a_hf1_score": [0.5]}
    )
    art = ucsc_track.build_track(
        gdf, prefix, str(tmp_path), run_tools=False, enzyme="encas12a", genome_fa=fa, encas12a_list_floor=0.0
    )
    row = open(art["bed"]).readline().rstrip("\n").split("\t")

    # pilot formula, recomputed here from the same inputs
    oc = 1
    c06 = sum(cfd for _c, _p, _s, cfd, pam in ots if weights[pam][1])
    w06 = sum(cfd * weights[pam][0] for _c, _p, _s, cfd, pam in ots if weights[pam][1])
    c03 = sum(cfd for _c, _p, _s, cfd, pam in ots if weights[pam][2])
    w03 = sum(cfd * weights[pam][0] for _c, _p, _s, cfd, pam in ots if weights[pam][2])
    f = _as_index(art["autosql"])
    # SPEC_DECIMALS is 3, so the written value can sit half a place (5e-4) off the
    # exact formula. Tightening this below that tolerance tests the rounding, not
    # the formula.
    tol = 5e-4
    assert float(row[f["enCas12a_spec_fa06"]]) == pytest.approx(1 / (oc + c06), abs=tol)
    assert float(row[f["enCas12a_spec_fa06_wt"]]) == pytest.approx(1 / (oc + w06), abs=tol)
    assert float(row[f["enCas12a_spec_fa03"]]) == pytest.approx(1 / (oc + c03), abs=tol)
    assert float(row[f["enCas12a_spec_fa03_wt"]]) == pytest.approx(1 / (oc + w03), abs=tol)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))



def test_ascas12a_specolor_is_additive(tmp_path):
    """_specColor is an ADDITIONAL specificity-only coloring; itemRgb is unchanged.

    Cas9 keeps a blended default in reserved/itemRgb (grey by MIT specificity, else a
    Doench-2016 ramp) and offers two alternative colourings; this track's are
    colour the same guide. So adding _specColor must not repaint the default view.
    """
    prefix = str(tmp_path / "off.csv")
    _write_ascas12a_enum(prefix)
    art = ucsc_track.build_track(
        _ascas12a_guide_df(), prefix, str(tmp_path), run_tools=False, enzyme="ascas12a"
    )
    rows = [l.rstrip("\n").split("\t") for l in open(art["bed"])]
    aslines = open(art["autosql"]).read()

    # _specColor is declared immediately after itemRgb, before guideSeq (Cas9 order)
    assert (aslines.index("itemRgb") < aslines.index("_crisprDeepCpf1Color")
            < aslines.index("_specColor") < aslines.index("guideSeq"))

    # field 9 = itemRgb, field 10 = _specColor; both are R,G,B triples
    for r in rows:
        for idx in (8, 9, 10):
            assert re.fullmatch(r"\d{1,3},\d{1,3},\d{1,3}", r[idx]), r[idx]

    # itemRgb is untouched by the two new colour columns: it still greys for close
    # and non-unique off-targets exactly as before.
    for r in rows:
        if r[9] in ("150,150,150", "120,120,120"):
            assert r[8] == r[9]


def test_spec_color_helper():
    """_spec_color is a pure function of the specificity/off-target axis."""
    import numpy as np
    nan = float("nan")
    z = np.zeros(4, dtype=int)
    # Cas9 cuts on specificity x100: <=50 red, <=70 olive, >70 green
    spec = np.array([0.90, 0.71, 0.70, 0.51, 0.50, 0.20, nan])
    z7 = np.zeros(7, dtype=int)
    assert list(ucsc_track._spec_color(np.array([False]*7), z7, z7, z7, spec)) == [
        "0,255,0", "0,255,0", "128,128,0", "128,128,0", "255,0,0", "255,0,0", "0,255,0"]
    # non-unique takes the worst tier but is still COLOURED, as Cas9 does
    assert list(ucsc_track._spec_color(np.array([True]*7), z7, z7, z7, spec)) == ["255,0,0"]*7
    # a close (1-mm) off-target does NOT grey this column (specificity alone)
    one7 = np.ones(7, dtype=int)
    assert ucsc_track._spec_color(np.array([False]*7), z7, one7, z7, spec)[0] == "0,255,0"


def test_cas9_style_field_names():
    assert ucsc_track._cas9_style_name("deepcpf1_score") == "deepcpf1"
    assert ucsc_track._cas9_style_name("ascas12a_deepcpf1_score") == "ascas12aDeepcpf1"
    assert ucsc_track._cas9_style_name("enpam_gb") == "enpamGb"


def test_efficiency_ramp_matches_cas9_cutoffs():
    """The red/yellow/green cuts follow the Cas9 crisprAllTargets description:
    red <= 30, yellow > 30 and <= 55, green > 55 (percentile). These were 50/75,
    which painted a percentile-60 guide yellow where Cas9 paints it green."""
    import numpy as np
    spec = np.full(7, 1.0)  # fully specific, so the ramp is what shows
    eff = np.array([0.0, 30.0, 30.1, 55.0, 55.1, 90.0, float("nan")])
    assert list(ucsc_track._spec_eff_color(spec, eff)) == [
        "255,100,100",  # 0   -> red
        "255,100,100",  # 30  -> red (boundary inclusive)
        "255,255,0",    # 30.1-> yellow
        "255,255,0",    # 55  -> yellow (boundary inclusive)
        "0,200,0",      # 55.1-> green
        "0,200,0",      # 90  -> green
        "0,200,0",      # NaN -> treated as fully efficient, as before
    ]
    assert (ucsc_track.EFF_RAMP_LOW, ucsc_track.EFF_RAMP_HIGH) == (30, 55)


# ---------------------------------------------------------------------------
# genome-wide percentile re-grading (chunked merges)
# ---------------------------------------------------------------------------
def test_rank_pct_matches_pct():
    """_rank_pct over a value histogram == _pct over the values themselves."""
    vals = [0.1, 0.25, 0.5, 0.5, 0.75, 1.0, 1.0, 1.0]
    keys = [f"{round(v, 4)}".encode() for v in vals]
    counts = {}
    for k in keys:
        counts[k] = counts.get(k, 0) + 1
    lut = ucsc_track._rank_pct(counts)
    ref = ucsc_track._pct(pd.Series(vals)).tolist()
    assert [lut[k] for k in keys] == pytest.approx(ref)


def test_regrade_bed_percentiles_is_genome_wide(tmp_path):
    """Two chunks whose score ranges barely overlap: per-chunk ranking makes the
    top of chunk A and the top of chunk B both look like the 75th percentile.
    After re-grading, every percentile is a rank over both chunks together, the
    ramp colours follow, greys survive and the hover is rewritten."""
    fields = [
        "chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd",
        "itemRgb", "_crisprDeepCpf1Color", "_specColor", "guideSeq", "pam", "unique",
        "AsCas12a_specificity", "deepcpf1", "deepcpf1Perc", "ascas12aDeepcpf1",
        "ascas12aDeepcpf1Perc", "_mouseOver", "_offset",
    ]
    spec = ucsc_track.track_regrade_spec("ascas12a", fields)

    def row(chrom, pos, spec_raw, dp, dp_pct, asc, asc_pct, rgb, alt, uniq="True"):
        return "\t".join([
            chrom, str(pos), str(pos + 27), "", "100", "+", str(pos + 4), str(pos + 27),
            rgb, alt, "0,255,0", "A" * 23, "TTTA", uniq, spec_raw,
            dp, dp_pct, asc, asc_pct, "stale hover", "0",
        ])

    # chunk A holds the low scores, chunk B the high ones; each ranked 25/50/75
    # within itself, so "75" means two different things.
    rows = [
        row("chr1", 100, "1.0", "10", "25", "0.10", "25", "255,100,100", "255,100,100"),
        row("chr1", 200, "1.0", "20", "50", "0.20", "50", "255,100,100", "255,100,100"),
        row("chr1", 300, "1.0", "30", "75", "0.30", "75", "255,100,100", "255,100,100"),
        row("chr2", 100, "1.0", "70", "25", "0.70", "25", "255,100,100", "255,100,100"),
        row("chr2", 200, "1.0", "80", "50", "0.80", "50", "255,100,100", "255,100,100"),
        row("chr2", 300, "1.0", "90", "75", "0.90", "75", "255,100,100", "255,100,100"),
        # a greyed, non-unique guide: the grey must survive and the hover must say so
        row("chr3", 100, "", "95", "75", "0.95", "75", "150,150,150", "150,150,150", uniq="False"),
    ]
    bed_in = tmp_path / "in.bed"
    bed_in.write_text("\n".join(rows) + "\n")
    bed_out = tmp_path / "out.bed"

    counts = {f: {} for f in [r for r, _, _, _ in spec["models"]]}
    for r in rows:
        f = r.split("\t")
        for name, i in (("deepcpf1", 15), ("ascas12aDeepcpf1", 17)):
            v = f[i].encode()
            if v:
                counts[name][v] = counts[name].get(v, 0) + 1

    n = ucsc_track._regrade_bed_percentiles(str(bed_in), str(bed_out), fields, spec, counts)
    assert n == len(rows)
    out = [l.split("\t") for l in bed_out.read_text().splitlines()]

    # 7 distinct scores, so the k-th smallest (1-indexed) ranks at 100*k/7
    expect = [round(100.0 * k / 7) for k in range(1, 8)]
    assert [int(r[16]) for r in out] == expect
    assert [int(r[18]) for r in out] == expect
    # the two "75"s have separated: chunk A's top is now mid-pack, chunk B's is top
    assert int(out[2][16]) < 50 < int(out[5][16])

    # ramp colours follow the new percentiles (<=30 red, <=55 yellow, >55 green)
    def ramp(pct):
        return "255,100,100" if pct <= 30 else ("255,255,0" if pct <= 55 else "0,200,0")

    assert [r[8] for r in out[:6]] == [ramp(p) for p in expect[:6]]
    assert [r[9] for r in out[:6]] == [ramp(p) for p in expect[:6]]
    # the grey guide keeps BOTH greys even though its raw score is the highest
    assert out[6][8] == "150,150,150" and out[6][9] == "150,150,150"
    # _specColor is an absolute cutoff, not a rank, so it is never touched
    assert {r[10] for r in out} == {"0,255,0"}

    # hover is regenerated, not left stale
    assert all("stale" not in r[19] for r in out)
    # the hover specificity is the raw score column read straight back, so it agrees
    # with column 5 (100 here) instead of reporting the mid-rank of the 1.0 tie
    assert out[5][4] == "100"
    assert out[5][19] == f"AsCas12a Spec. Score: 100, Kim 2018: {expect[5]}%, Chen 2025: {expect[5]}%"
    assert out[6][19] == "Not unique, Kim 2018: 100%, Chen 2025: 100%"


@pytest.mark.parametrize(
    "decimals,values,expect",
    [
        # seq-DeepCpf1 shape, 0 to 1 at 2 decimals: 1.0 and 0.99996 both show as 1.00
        (2, [1.0, 0.99996, 0.994, 0.5, 0.5, 0.123456, float("nan")],
         ["1.00", "1.00", "0.99", "0.50", "0.50", "0.12", ""]),
        # DeepCpf1 shape, about 1 to 115 as an integer: 50.9896 and 51.0023 both show as 51
        (0, [50.9896, 51.0023, 114.672, 0.9584, 0.9584, 30.4, float("nan")],
         ["51", "51", "115", "1", "1", "30", ""]),
    ],
)
def test_percentiles_rank_the_displayed_value_in_both_paths(decimals, values, expect):
    """build_track and the merge re-grade must agree exactly, at either precision.

    build_track ranks in pandas over a float column; the merge re-grade can only see
    the strings written to the BED. Ranking the full-precision float in one and the
    written value in the other makes a chunked build disagree with an unchunked one
    wherever rounding merges distinct floats, so both rank the DISPLAYED value.

    That is also why reducing a model's declared decimals re-ranks it rather than
    merely reformatting it: the two values that tie here are separated by a full-float
    rank, and the last assertion pins exactly that."""
    raw = pd.Series(values)
    disp = ucsc_track._fmt_score(raw, decimals)
    assert list(disp) == expect

    from_build = ucsc_track._pct_of_displayed(raw, decimals).tolist()
    counts = {}
    for v in disp:
        if v:
            counts[v.encode()] = counts.get(v.encode(), 0) + 1
    lut = ucsc_track._rank_pct(counts)
    from_merge = [None if not v else lut[v.encode()] for v in disp]

    assert from_build[:-1] == pytest.approx(from_merge[:-1])
    assert pd.isna(from_build[-1]) and from_merge[-1] is None
    # the two values that round together tie, where a full-float rank separates them
    assert from_build[0] == from_build[1]
    assert ucsc_track._pct(raw).tolist()[0] != ucsc_track._pct(raw).tolist()[1]


def test_hover_specificity_is_raw_not_a_percentile(tmp_path):
    """The hover specificity must equal column 5 (raw x100), as on the Cas9 track,
    whose hover reads "MIT Spec. Score: 80" with score == 80.

    A rank is wrong here: specificity is already absolute, and its distribution is
    dominated by the perfectly-specific tie, so ranking it reports a low number for
    a guide the score column calls 100."""
    prefix = str(tmp_path / "off.csv")
    _write_enum(prefix)
    art = ucsc_track.build_track(
        _guide_df(), prefix, str(tmp_path), run_tools=False, enzyme="ascas12a"
    )
    rows = [l.split("\t") for l in open(art["bed"]).read().splitlines()]
    assert rows
    for r in rows:
        # _mouseOver is always the second-to-last field, _offset the last, so this
        # does not depend on how many on-target models the fixture carries
        score, hover = r[4], r[-2]
        if hover.startswith("Not unique"):
            continue
        assert hover.startswith(f"AsCas12a Spec. Score: {score},") or hover == f"AsCas12a Spec. Score: {score}", (
            f"hover {hover!r} disagrees with score {score!r}"
        )


def test_track_shape_is_universal_across_enzymes():
    """Every enzyme ships the same track shape. The only legitimate difference is
    which specificity columns and which cleavage models it carries, so the colour
    block, the split raw/...Perc layout and the trailer must be identical."""
    import re

    def names(text):
        return [m.group(1) for m in re.finditer(
            r"^\s+(?:string|uint|char|bigint)\S*\s+(\w+);", text, re.M)]

    a = names(ucsc_track.build_autosql_ascas12a(
        ["deepcpf1_score", "ascas12a_deepcpf1_score"]))
    e = names(ucsc_track.build_autosql_encas12a(
        ["enpam_gb_score", "encas12a_hf1_score"]))

    # bed9 + the three colourings, in the same order, before any enzyme-specific column
    bed9 = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
            "thickStart", "thickEnd"]
    head = bed9 + list(ucsc_track.COLOR_FIELDS) + ["guideSeq", "pam", "unique"]
    assert a[:len(head)] == head
    assert e[:len(head)] == head

    # same trailer
    assert a[-2:] == ["_mouseOver", "_offset"]
    assert e[-2:] == ["_mouseOver", "_offset"]

    # each cleavage model is a raw column plus a ...Perc sibling, never a packed string
    for fields, enzyme in ((a, "ascas12a"), (e, "encas12a")):
        spec = ucsc_track.track_regrade_spec(enzyme, fields)
        assert spec is not None, enzyme
        assert len(spec["models"]) == 2
        for raw, _label, colour, _dec in spec["models"]:
            assert raw in fields
            assert raw + "Perc" in fields
            assert colour in ucsc_track.EFF_COLOR_FIELDS
        # first model drives itemRgb, second the alternative colour
        assert [c for _, _, c, _ in spec["models"]] == list(ucsc_track.EFF_COLOR_FIELDS)
        assert spec["spec_score"] == "score"


def test_regrade_spec_follows_the_models_actually_scored():
    """An extra or missing cleavage score changes only how many models are graded,
    never how any of them is graded."""
    one = ucsc_track.build_autosql_encas12a(["enpam_gb_score"])
    import re
    f = [m.group(1) for m in re.finditer(
        r"^\s+(?:string|uint|char|bigint)\S*\s+(\w+);", one, re.M)]
    spec = ucsc_track.track_regrade_spec("encas12a", f)
    assert [raw for raw, _, _, _ in spec["models"]] == ["enpamGb"]
    assert [c for _, _, c, _ in spec["models"]] == ["itemRgb"]
    assert "encas12aHf1" not in f

    # an enzyme with no entry has no rule; build_track rejects it outright
    assert ucsc_track.track_regrade_spec("nosuchenzyme", f) is None


def test_marginal_specificity_is_not_greyed():
    """A 0.20-to-0.50 specificity guide shows its efficiency ramp, not a grey.

    The old scheme had a second grey tier there, which hid the efficiency of a guide
    whose specificity is poor but usable. Only below SPEC_GREY_BAD is a unique,
    no-close-off-target guide greyed now.
    """
    import numpy as np

    spec = np.array([0.10, 0.19, 0.20, 0.35, 0.49, 0.50, 0.95])
    top = np.full(spec.shape, 90.0)      # high efficiency -> green if not greyed
    bottom = np.full(spec.shape, 5.0)    # low efficiency  -> red if not greyed

    assert list(ucsc_track._spec_eff_color(spec, top)) == [
        "50,50,50", "50,50,50", "0,200,0", "0,200,0", "0,200,0", "0,200,0", "0,200,0",
    ]
    # the ramp, not a grey, is what varies across the 0.20-0.50 band
    assert list(ucsc_track._spec_eff_color(spec, bottom)) == [
        "50,50,50", "50,50,50", "255,100,100", "255,100,100", "255,100,100",
        "255,100,100", "255,100,100",
    ]
    assert "90,90,90" not in set(ucsc_track._spec_eff_color(spec, top))
    assert not hasattr(ucsc_track, "SPEC_GREY_MARGINAL")


def test_regrade_converts_the_retired_marginal_grey(tmp_path):
    """Re-grading an older track repaints its 90,90,90 cells with the ramp, while the
    uniqueness and bad-specificity greys survive untouched.

    This is what lets a finished track adopt the new scheme with no rescan: the colour
    is a function of the percentile, which the track already carries.
    """
    fields = [
        "chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart",
        "thickEnd", "itemRgb", "_crisprDeepCpf1Color", "_specColor", "guideSeq", "pam",
        "unique", "AsCas12a_specificity", "deepcpf1", "deepcpf1Perc", "ascas12aDeepcpf1",
        "ascas12aDeepcpf1Perc", "_mouseOver", "_offset",
    ]
    spec = ucsc_track.track_regrade_spec("ascas12a", fields)

    def row(dp, rgb, alt, uniq="True"):
        return "\t".join([
            "chrI", "100", "127", "", "30", "+", "104", "127", rgb, alt, "255,0,0",
            "A" * 23, "TTTA", uniq, "0.3000", dp, "0", dp, "0", "stale", "0",
        ])

    src = tmp_path / "in.bed"
    src.write_text("\n".join([
        row("90.0", "90,90,90", "90,90,90"),        # retired marginal grey, high eff
        row("10.0", "90,90,90", "90,90,90"),        # retired marginal grey, low eff
        row("90.0", "50,50,50", "50,50,50"),        # bad specificity: stays grey
        row("90.0", "120,120,120", "120,120,120"),  # close off-target: stays grey
        row("90.0", "150,150,150", "150,150,150", "False"),  # non-unique: stays grey
    ]) + "\n")

    counts = {"deepcpf1": {b"90.0": 4, b"10.0": 1},
              "ascas12aDeepcpf1": {b"90.0": 4, b"10.0": 1}}
    out = tmp_path / "out.bed"
    n = ucsc_track._regrade_bed_percentiles(str(src), str(out), fields, spec, counts)
    assert n == 5

    got = [l.split("\t") for l in out.read_text().splitlines()]
    # the two marginal-grey rows now carry a ramp colour driven by their percentile
    assert got[0][8] == "0,200,0" and got[0][9] == "0,200,0"
    assert got[1][8] == "255,100,100" and got[1][9] == "255,100,100"
    # the three real greys are preserved
    assert got[2][8] == "50,50,50"
    assert got[3][8] == "120,120,120"
    assert got[4][8] == "150,150,150"
    # nothing repainted a cell to the retired grey
    assert all(r[8] != "90,90,90" and r[9] != "90,90,90" for r in got)


ASCAS12A_FIELDS = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart",
    "thickEnd", "itemRgb", "_crisprDeepCpf1Color", "_specColor", "guideSeq", "pam",
    "unique", "AsCas12a_specificity", "deepcpf1", "deepcpf1Perc", "ascas12aDeepcpf1",
    "ascas12aDeepcpf1Perc", "_mouseOver", "_offset",
]


def _valid_row(**over):
    """One row that satisfies every invariant; keyword args break exactly one."""
    f = {
        "chrom": "chrI", "chromStart": "100", "chromEnd": "127", "name": "",
        "score": "90", "strand": "+", "thickStart": "104", "thickEnd": "127",
        "itemRgb": "0,200,0", "_crisprDeepCpf1Color": "0,200,0", "_specColor": "0,255,0",
        "guideSeq": "A" * 23, "pam": "TTTA", "unique": "True",
        # each cleavage score at the precision its model declares in
        # ENZYME_TRACK_SPECS: DeepCpf1 as an integer, the 0-to-1 model at 2 decimals
        "AsCas12a_specificity": "0.9000", "deepcpf1": "70", "deepcpf1Perc": "80",
        "ascas12aDeepcpf1": "0.70", "ascas12aDeepcpf1Perc": "80",
        "_mouseOver": "AsCas12a Spec. Score: 90, Kim 2018: 80%, Chen 2025: 80%",
        "_offset": "0",
    }
    f.update(over)
    return "\t".join(f[k] for k in ASCAS12A_FIELDS)


def _write_bed(tmp_path, rows, name="t.bed"):
    p = tmp_path / name
    p.write_text("\n".join(rows) + "\n")
    return str(p)


def test_validator_accepts_a_well_formed_track(tmp_path):
    bed = _write_bed(tmp_path, [_valid_row(), _valid_row(chromStart="200", chromEnd="227")])
    stats = ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")
    assert stats["guides"] == 2


def test_validator_catches_chunk_local_percentiles(tmp_path):
    """The regression that shipped: one raw score mapping to two percentiles, which is
    exactly the signature of ranking within a chunk instead of genome-wide."""
    bed = _write_bed(tmp_path, [
        _valid_row(deepcpf1="70", deepcpf1Perc="80"),
        _valid_row(chromStart="200", chromEnd="227", deepcpf1="70", deepcpf1Perc="42"),
    ])
    with pytest.raises(ValueError, match="not genome-wide"):
        ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")


def test_validator_catches_a_retired_or_stray_colour(tmp_path):
    """The other regression that shipped: a colour the scheme no longer produces."""
    bed = _write_bed(tmp_path, [_valid_row(itemRgb="90,90,90")])
    with pytest.raises(ValueError, match="cannot produce"):
        ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")


def test_validator_catches_score_disagreeing_with_specificity(tmp_path):
    bed = _write_bed(tmp_path, [_valid_row(score="42")])
    with pytest.raises(ValueError, match="differ by more than 1"):
        ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")


def test_validator_tolerates_the_four_decimal_rounding_of_specificity(tmp_path):
    """score is rounded from the full-precision float, the column is written to 4
    decimals, so a value sitting exactly on a half legitimately differs by 1. Real on
    sacCer3: 113 of 590,425 guides, e.g. a displayed 0.905 whose score is 91."""
    ok = _write_bed(tmp_path, [
        _valid_row(AsCas12a_specificity="0.9050", score="91",
                   _mouseOver="AsCas12a Spec. Score: 91, Kim 2018: 80%, Chen 2025: 80%"),
        _valid_row(chromStart="200", chromEnd="227",
                   AsCas12a_specificity="0.9150", score="91",
                   _mouseOver="AsCas12a Spec. Score: 91, Kim 2018: 80%, Chen 2025: 80%"),
    ], name="ok.bed")
    assert ucsc_track.validate_track_bed(ok, ASCAS12A_FIELDS, "ascas12a")["guides"] == 2

    # two off is drift, not rounding
    bad = _write_bed(tmp_path, [
        _valid_row(AsCas12a_specificity="0.9000", score="92",
                   _mouseOver="AsCas12a Spec. Score: 92, Kim 2018: 80%, Chen 2025: 80%")
    ], name="bad.bed")
    with pytest.raises(ValueError, match="differ by more than 1"):
        ucsc_track.validate_track_bed(bad, ASCAS12A_FIELDS, "ascas12a")


def test_validator_catches_a_percentile_in_the_hover(tmp_path):
    """The hover must repeat `score`, never a rank of the specificity."""
    bed = _write_bed(tmp_path, [
        _valid_row(_mouseOver="AsCas12a Spec. Score: 53, Kim 2018: 80%, Chen 2025: 80%")
    ])
    with pytest.raises(ValueError, match="hover specificity"):
        ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")


def test_validator_catches_a_ramp_colour_on_a_non_unique_guide(tmp_path):
    bed = _write_bed(tmp_path, [_valid_row(unique="False", itemRgb="0,200,0")])
    with pytest.raises(ValueError, match="ramp colour instead of a grey"):
        ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")


def test_validator_catches_a_percentile_without_a_raw_value(tmp_path):
    bed = _write_bed(tmp_path, [_valid_row(deepcpf1="", deepcpf1Perc="80")])
    with pytest.raises(ValueError, match="without a percentile, or the reverse"):
        ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")


def test_validator_catches_a_wrong_field_count(tmp_path):
    bed = _write_bed(tmp_path, [_valid_row() + "\textra"])
    with pytest.raises(ValueError, match="autoSql declares"):
        ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")


def test_validator_catches_a_score_at_the_wrong_precision(tmp_path):
    """A cleavage score written to more places than its model declares fails the build.

    This is the check that stops a display-only "rounding" change from shipping: a
    column written at one precision and ranked at another breaks the raw-to-percentile
    invariant, and the two checks together make that unshippable rather than something
    to notice later on the browser.
    """
    bed = _write_bed(tmp_path, [_valid_row(deepcpf1="70.0")])
    with pytest.raises(ValueError, match="not written to the 0 decimals"):
        ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")

    bed = _write_bed(tmp_path, [_valid_row(ascas12aDeepcpf1="0.7345")], name="b.bed")
    with pytest.raises(ValueError, match="not written to the 2 decimals"):
        ucsc_track.validate_track_bed(bed, ASCAS12A_FIELDS, "ascas12a")


def test_declared_decimals_are_the_model_scale_not_a_house_style():
    """Each model's precision has to suit its own output range.

    DeepCpf1 spans roughly 1 to 115, so an integer still separates guides; the three
    0-to-1 models would collapse to 2 distinct values as integers and need 2 places.
    Measured distinct-value counts: DeepCpf1 112 as integers, and seq-DeepCpf1 /
    EnPAM-GB / enAsCas12a-HF1 92 to 94 at 2 decimals, all inside the 92 to 126 band
    the Cas9 track's six score fields occupy. Pinned because a uniform precision
    across models is the tempting wrong answer.
    """
    dec = {
        col: d
        for espec in ucsc_track.ENZYME_TRACK_SPECS.values()
        for col, _lab, d in espec["ontarget"]
    }
    assert dec["deepcpf1_score"] == 0
    assert dec["ascas12a_deepcpf1_score"] == 2
    assert dec["enpam_gb_score"] == 2
    assert dec["encas12a_hf1_score"] == 2


def test_rounding_a_score_re_ranks_it_rather_than_only_reformatting_it():
    """Coarsening the display changes the percentile, which is the whole reason this
    is a track change and not a cosmetic one.

    Three floats that differ at the fourth decimal rank apart on the full value and
    tie once displayed as integers, so they must all carry the same percentile."""
    raw = pd.Series([50.9896, 51.0023, 51.4, 10.0])
    assert list(ucsc_track._fmt_score(raw, 0)) == ["51", "51", "51", "10"]

    pct = ucsc_track._pct_of_displayed(raw, 0).tolist()
    assert pct[0] == pct[1] == pct[2], "values that display alike must rank alike"
    assert pct[3] < pct[0]
    # at 4 decimals they would have been three different percentiles
    assert len(set(ucsc_track._pct_of_displayed(raw, 4).tolist())) == 4


def test_specificity_is_three_decimals_keeping_one_trailing_zero():
    """Specificity displays at 3 places, trailing zeros stripped but one kept.

    A CFD aggregate does not support a fourth place and the Cas9 track shows none.
    Rounding it is safe where rounding a cleavage score is not, because specificity is
    never ranked, so no percentile or colour moves with it. An exact 1 must still read
    `1.0`, not `1` and not `1.000`.
    """
    assert ucsc_track.SPEC_DECIMALS == 3
    raw = pd.Series([1.0, 0.1601, 0.9050, 0.1, 0.0, 0.12345, float("nan")])
    assert list(ucsc_track._fmt_raw(raw)) == ["1.0", "0.16", "0.905", "0.1", "0.0", "0.123", ""]

    # worst case is half of the last place, i.e. 0.05 of a `score` unit, which has to
    # stay well inside the tolerance validate_track_bed allows between the two
    import numpy as np
    x = np.linspace(0, 1, 10001)
    disp = ucsc_track._fmt_raw(pd.Series(x)).astype(float).to_numpy()
    assert np.max(np.abs(disp - x)) <= 0.0005 + 1e-12
    assert np.max(np.abs(disp - x)) * 100.0 < 1.0


def test_hover_follows_the_cas9_shape_and_cites_the_right_model(tmp_path):
    """The hover must read like Cas9's "MIT Spec. Score: 48, Doench 2016: 18%,
    Moreno-Mateos: 72%": a "Spec. Score" label, then one cited model per cleavage
    figure, each a percentile with a % suffix.

    "DeepCpf1" is specifically wrong as a label. crisprware loads
    weights/Seq_deepCpf1_weights.h5, the sequence-only model; DeepCpf1 proper is the
    chromatin-accessibility variant and is not what runs. Tracks shipped before
    2026-10-19 carried the wrong name, so this pins the citation form.
    """
    prefix = str(tmp_path / "off.csv")
    _write_enum(prefix)
    art = ucsc_track.build_track(
        _guide_df(), prefix, str(tmp_path), run_tools=False, enzyme="ascas12a"
    )
    rows = [l.split("\t") for l in open(art["bed"]).read().splitlines()]
    assert rows
    import re
    for r in rows:
        hover = r[-2]
        if hover.startswith("Not unique"):
            continue
        assert hover.startswith("AsCas12a Spec. Score: "), hover
        assert "DeepCpf1" not in hover, "the model that runs is Seq-deepCpf1, not DeepCpf1"
        for label in ("Kim 2018", "Chen 2025"):
            if label + ":" in hover:
                assert re.search(rf"{label}: \d+%", hover), hover

    labels = [lab for _c, lab, _d in ucsc_track.ENZYME_TRACK_SPECS["ascas12a"]["ontarget"]]
    assert labels == ["Kim 2018", "Chen 2025"]
    assert ucsc_track.ENZYME_TRACK_SPECS["ascas12a"]["spec_label"] == "AsCas12a Spec. Score"
