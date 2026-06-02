"""Unit tests for crisprware.ucsc_track (UCSC Cas12a track assembly).

Exercises the adapter on synthetic crispr-ots streaming output: the strand-aware
5'-PAM bigBed geometry, the off-target browser-pos mapping (+: pos=start,
-: pos=start+3), _mismatchCounts/_crisprOfftargets assembly, and the _offset ->
byte-offset round-trip. No TensorFlow / GPU / UCSC tools required (run_tools=False).
"""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from crisprware import ucsc_track  # noqa: E402


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
    art = ucsc_track.build_track(_guide_df(), prefix, str(tmp_path), run_tools=False)

    bed = [line.rstrip("\n").split("\t") for line in open(art["bed"])]
    assert len(bed) == 2
    by_chromstart = {int(r[1]): r for r in bed}

    # g_uniq (+): spacer [104,127) -> site [100,127) (27nt), thick [104,127) (23nt)
    plus = by_chromstart[100]
    assert int(plus[2]) - int(plus[1]) == 27  # chromEnd - chromStart
    assert int(plus[7]) - int(plus[6]) == 23  # thickEnd - thickStart
    assert plus[6] == "104" and plus[7] == "127"  # thick = spacer
    assert plus[10] == "TTTA"  # 4-nt PAM
    assert plus[5] == "+"
    assert int(plus[-1]) > 0  # _offset set (has off-targets)

    # g_dup (-): spacer [204,227) -> site [204,231) (27nt), thick [204,227) (PAM on the right)
    minus = next(r for r in bed if r[5] == "-")
    assert int(minus[1]) == 204 and int(minus[2]) == 231  # chromStart, chromEnd
    assert int(minus[6]) == 204 and int(minus[7]) == 227  # thick = spacer
    assert minus[8] == "150,150,150"  # grey (non-unique/dropped)
    assert minus[-1] == "0"  # dropped -> no crisprDetails line

    # All geometry invariants
    for r in bed:
        assert int(r[2]) - int(r[1]) == 27
        assert int(r[7]) - int(r[6]) == 23
        assert len(r[10]) == 4


def test_offtarget_browser_pos_and_offset_seek(tmp_path):
    prefix = str(tmp_path / "off.csv")
    _write_enum(prefix)
    art = ucsc_track.build_track(_guide_df(), prefix, str(tmp_path), run_tools=False)

    # find g_uniq's _offset
    bed = [line.rstrip("\n").split("\t") for line in open(art["bed"])]
    off = next(int(r[-1]) for r in bed if r[1] == "100")

    raw = open(art["details_tab"], "rb").read()
    # header at byte 0; the recorded offset points at g_uniq's line
    line = raw[off:].split(b"\n", 1)[0].decode()
    counts, lst = line.split("\t")
    assert counts == "0,1,0,2,0"  # mismatch_counts ; -> ,
    # off-targets sorted by score desc: + at 500 (pos=start=500), - at 800 (pos=start+3=803)
    entries = lst.split("|")
    assert entries[0] == "chr1;500+;800"  # cfd 0.80 -> 800; + strand pos=start
    assert entries[1] == "chr3;803-;300"  # cfd 0.30 -> 300; - strand pos=start+3
    # sum of listed strong off-targets <= sum of mismatch counts (=3)
    assert len(entries) <= sum(int(x) for x in counts.split(","))


def test_blank_list_when_too_many(tmp_path):
    prefix = str(tmp_path / "off.csv")
    _write_enum(prefix)
    # off_target_count 3 > blank_threshold 1 -> list blanked, counts kept
    art = ucsc_track.build_track(_guide_df(), prefix, str(tmp_path), blank_threshold=1, run_tools=False)
    raw = open(art["details_tab"], "rb").read().decode()
    body = raw.split("\n", 1)[1]
    counts, lst = body.split("\n")[0].split("\t")
    assert counts == "0,1,0,2,0"
    assert lst == ""  # blanked


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
