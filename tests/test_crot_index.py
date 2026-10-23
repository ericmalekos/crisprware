"""The off-target index's enumerated PAM set, read off disk and checked per enzyme.

This exists because of a specific failure. AsCas12a does not cleave a TTTT-prefixed
protospacer, so a TTTT site in the off-target index is a false off-target. The scanner
excludes TTTT from exactly ONE of its outputs, ``specificity_tttv``; the mismatch
buckets behind the close-off-target grey, the off-target list in the details popup, the
off-target count, max CFD, the saturation cap and the 0-mm screen behind ``unique`` all
count every site they are handed. Three separate comments in the track builder asserted
the index was "TTTV-only" while eight genomes had in fact been indexed with TTTN, and
nothing in the stack could tell.

So the assumption is read out of the index header instead of asserted in prose, and
these tests pin both halves: the header parse (a wrong struct format returns an EMPTY
PAM list, which would silently pass a check whose job is to reject TTTT) and the
per-enzyme verdict.
"""
import struct

import pandas as pd
import pytest

from crisprware.utils import crot_index as ci


def _write_crot(path, pams, *, enzyme_name="PAM-TTTN-5p-23nt", magic=ci.MAGIC,
                proto_len=23, bin_width=14, pam_side=ci.PAM_FIVE_PRIME, entries=945116):
    """Write a synthetic .crot carrying only the header and its string sections.

    The layout mirrors ``DbHeader`` in crispr-db/src/mmap_db.rs: a 128-byte repr(C)
    little-endian header, then the enzyme name, then ``pam_count * pam_len`` ASCII
    bytes of expanded PAMs. Writing it here rather than shelling out to
    ``crispr-ots index`` keeps the test hermetic AND documents the format, so a
    change on the Rust side shows up as a failure rather than as an empty PAM list.
    """
    name = enzyme_name.encode()
    pam_len = len(pams[0]) if pams else 4
    blob = "".join(pams).encode()
    enzyme_name_off = ci.HEADER_BYTES
    pam_data_off = enzyme_name_off + len(name)
    contig_off = pam_data_off + len(blob)
    header = struct.pack(
        "<8s10I10Q",
        magic, 2, 0, bin_width, pam_side, pam_len, proto_len,
        len(pams), 0, 0, 0,
        entries, 0, 0, enzyme_name_off, pam_data_off, contig_off, contig_off, contig_off, 0, 0,
    )
    assert len(header) == ci.HEADER_BYTES
    path.write_bytes(header + name + blob)
    return path


# --- the parse ---------------------------------------------------------------

def test_header_round_trips_every_field(tmp_path):
    p = _write_crot(tmp_path / "i.crot", ["TTTA", "TTTC", "TTTG", "TTTT"])
    h = ci.read_crot_header(p)
    assert h.pams == ("TTTA", "TTTC", "TTTG", "TTTT")
    assert h.enzyme_name == "PAM-TTTN-5p-23nt"
    assert (h.protospacer_len, h.bin_width, h.total_entries) == (23, 14, 945116)
    assert h.pam_side_label == "5'"


def test_struct_format_is_exactly_the_header_size():
    # The load-bearing invariant. Get this wrong and enzyme_name_off / pam_data_off
    # read as 0, read_crot_header returns pams == (), and has_tttt() is False for a
    # TTTN index, i.e. the check passes on exactly the input it exists to reject.
    assert ci._H.size == ci.HEADER_BYTES == 128


def test_a_stem_and_an_explicit_crot_resolve_to_the_same_file(tmp_path):
    p = _write_crot(tmp_path / "i.crot", ["TTTA"])
    assert ci.resolve_index_path(tmp_path / "i") == p
    assert ci.resolve_index_path(p) == p
    # The stem crispr-ots is given on the command line carries the directory name
    # twice (<stem>_crisprots/<stem>_crisprots), and has no suffix to replace.
    assert ci.resolve_index_path("/x/y_crisprots/y_crisprots").name == "y_crisprots.crot"


def test_non_crispr_ots_file_is_rejected(tmp_path):
    p = _write_crot(tmp_path / "bad.crot", ["TTTA"], magic=b"NOTACROT")
    with pytest.raises(ValueError, match="not a crispr-ots v2 index"):
        ci.read_crot_header(p)


def test_truncated_file_is_rejected(tmp_path):
    p = tmp_path / "short.crot"
    p.write_bytes(b"CRSPROT2" + b"\0" * 40)
    with pytest.raises(ValueError, match="shorter than a crispr-ots header"):
        ci.read_crot_header(p)


# --- the per-enzyme verdict --------------------------------------------------

def test_ascas12a_refuses_a_tttn_index(tmp_path):
    p = _write_crot(tmp_path / "tttn.crot", ["TTTA", "TTTC", "TTTG", "TTTT"])
    with pytest.raises(ValueError, match="enumerates TTTT"):
        ci.assert_enzyme_pams(p, "ascas12a")


def test_ascas12a_accepts_a_tttv_index(tmp_path):
    p = _write_crot(tmp_path / "tttv.crot", ["TTTA", "TTTC", "TTTG"],
                    enzyme_name="PAM-TTTV-5p-23nt", entries=590426)
    h = ci.assert_enzyme_pams(p, "ascas12a")
    assert not h.has_tttt()
    assert h.total_entries == 590426


def test_encas12a_keeps_a_tttn_index(tmp_path):
    # enAsCas12a DOES cleave TTTT, so excluding it would be the error here.
    p = _write_crot(tmp_path / "tttn.crot", ["TTTA", "TTTC", "TTTG", "TTTT"])
    assert ci.assert_enzyme_pams(p, "encas12a").has_tttt()


def test_the_refusal_names_the_rebuild_command(tmp_path):
    # The message has to be actionable from a shard log, where the reader has no
    # context beyond one failing task.
    p = _write_crot(tmp_path / "tttn.crot", ["TTTA", "TTTC", "TTTG", "TTTT"])
    with pytest.raises(ValueError) as e:
        ci.assert_enzyme_pams(p, "ascas12a")
    msg = str(e.value)
    assert "--pam TTTV" in msg and "-l 23" in msg and "--bin_width 14" in msg


def test_a_lowercase_pam_is_still_tttt(tmp_path):
    p = _write_crot(tmp_path / "lc.crot", ["ttta", "tttt"])
    with pytest.raises(ValueError, match="enumerates TTTT"):
        ci.assert_enzyme_pams(p, "ascas12a")


# --- wired into the real entry point ----------------------------------------

def test_ucscgb_track_refuses_a_tttn_index_before_scanning(tmp_path):
    """The check that was missing: the track entry point itself must refuse.

    A unit-tested helper nobody calls is how this defect shipped, so this test goes
    through ``run_ucscgb_track`` and asserts it raises before any scan.
    """
    import argparse

    from crisprware import score_guides as sg

    idx = _write_crot(tmp_path / "g_crisprots.crot", ["TTTA", "TTTC", "TTTG", "TTTT"])
    outdir = tmp_path / "track"
    composite = "id,sequence,pam,chromosome,position,sense"
    gRNADF = pd.DataFrame({composite: ["chr1:1:+\tACGT\tTTTA\tchr1\t1\t+"]})
    args = argparse.Namespace(
        guidescan2_indices=[str(idx.with_suffix(""))],
        ucscgb=str(outdir),
        ucscgb_enzyme="ascas12a",
        ucscgb_reuse_offtargets=False,
    )
    with pytest.raises(ValueError, match="enumerates TTTT"):
        sg.run_ucscgb_track(args, gRNADF)
