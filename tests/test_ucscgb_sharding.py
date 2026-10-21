"""Phased (sharded) --ucscgb_chunk_max builds.

Chunk work was already independent: each chunk is a separate process with its own
output dir and no shared state. What was missing is a way to let something other than
one process run them, because the worker count is ``max(1, len(gpus))``, which is 1
on the CPU scanner. A 48-chunk genome then runs strictly serially.

These tests pin the two things that make a sharded build safe rather than fast:
the split is reproducible, and the merge refuses an incomplete set.
"""
import argparse
import os

import pytest

from crisprware import score_guides as sg


def _args(bed, outdir, chunk_max=10, **over):
    a = argparse.Namespace(
        grna_bed=str(bed),
        ucscgb=str(outdir),
        ucscgb_chunk_max=chunk_max,
        ucscgb_gpus=None,
        ucscgb_scanner="cpu",
        ucscgb_keep_chunks=True,
        ucscgb_phase="all",
        ucscgb_shard=None,
        ucscgb_enzyme="ascas12a",
        chrom_sizes="/dev/null",
        threads=1,
        output_directory=str(outdir),
    )
    for k, v in over.items():
        setattr(a, k, v)
    return a


def _write_bed(path, n):
    with open(path, "w") as fh:
        fh.write("#chrom\tstart\tstop\n")
        for i in range(n):
            fh.write(f"chr1\t{i * 100}\t{i * 100 + 27}\n")
    return path


def test_split_is_reproducible_and_keeps_every_guide(tmp_path):
    """The split must be a pure function of (guide BED, chunk_max).

    The merge job runs in a different process, possibly days later, and recomputes the
    expected chunk count from the BED rather than trusting what it finds on disk. That
    only works if the split is reproducible."""
    bed = _write_bed(tmp_path / "g.bed", 95)
    root = tmp_path / "chunks"
    root.mkdir()

    first = sg._ucscgb_split(_args(bed, tmp_path), str(root))
    names = [os.path.basename(p) for p in first]
    bodies = {p: open(p).read() for p in first}

    second = sg._ucscgb_split(_args(bed, tmp_path), str(root))
    assert [os.path.basename(p) for p in second] == names
    assert {p: open(p).read() for p in second} == bodies

    assert sg._ucscgb_expected_chunks(_args(bed, tmp_path)) == len(first)
    # every guide appears exactly once, and each chunk carries the header
    total = 0
    for p in first:
        lines = open(p).read().splitlines()
        assert lines[0].startswith("#chrom")
        total += len(lines) - 1
    assert total == 95


def test_chunk_paths_agree_across_phases(tmp_path):
    """A shard job and the merge job must derive the same paths for a chunk."""
    cb, tdir, odir = sg._ucscgb_chunk_paths(str(tmp_path), 7)
    assert cb.endswith("chunk_007.bed")
    assert tdir.endswith("track_007")
    assert odir.endswith("out_007")


def test_shard_refuses_without_a_split(tmp_path):
    """A shard must never split for itself: N concurrent array tasks would write the
    same chunk files at the same time and interleave their contents, giving chunks
    that are individually well-formed and collectively wrong."""
    bed = _write_bed(tmp_path / "g.bed", 50)
    a = _args(bed, tmp_path, ucscgb_phase="shard", ucscgb_shard=0)
    with pytest.raises(RuntimeError, match="run --ucscgb_phase split first"):
        sg.run_ucscgb_chunked(a)


def test_shard_index_must_be_in_range(tmp_path):
    bed = _write_bed(tmp_path / "g.bed", 50)  # chunk_max 10 -> 5 chunks
    a = _args(bed, tmp_path, ucscgb_phase="shard", ucscgb_shard=5)
    with pytest.raises(ValueError, match="out of range"):
        sg.run_ucscgb_chunked(a)
    a.ucscgb_shard = None
    with pytest.raises(ValueError, match="requires --ucscgb_shard"):
        sg.run_ucscgb_chunked(a)


# The basenames a real shard leaves in its chunk dir, read off
# ce11/track_ce11_cas12a_shardtest/chunks/track_000 after job 42472186_0 (2026-10-20).
# Hard-coding them here is the point: the first version of this fixture invented
# `ascas12a_crisprDetails.tab`, which matched the guard but not the build, so the test
# passed while the guard rejected four complete chunks.
REAL_CHUNK_FILES = (
    "ascas12a.bb",
    "ascas12a_crisprDetails.tab.gz",
    "ascas12a_crisprDetails.tab.gz.gzi",
    "ascas12aTargets.as",
    "guides.bed",
)


def _fake_chunk(root, i, as_text='table T\n"t"\n( string chrom; "c" )\n', files=REAL_CHUNK_FILES):
    _cb, tdir, _od = sg._ucscgb_chunk_paths(root, i)
    os.makedirs(tdir, exist_ok=True)
    for fn in files:
        open(os.path.join(tdir, fn), "w").write("x")
    open(os.path.join(tdir, "ascas12aTargets.as"), "w").write(as_text)
    return tdir


def test_merge_accepts_the_files_a_real_shard_leaves(tmp_path):
    """Every chunk of the ce11 shard test was complete and the merge refused all four
    (job 42472187, 2026-10-20): the guard required `ascas12a_crisprDetails.tab` while
    build_track writes it bgzipped. The required names have to be the ones
    merge_ucscgb_tracks opens, so a complete set must get past this check."""
    bed = _write_bed(tmp_path / "g.bed", 30)
    root = tmp_path / "chunks"
    root.mkdir()
    for i in range(3):
        _fake_chunk(str(root), i)

    a = _args(bed, tmp_path, ucscgb_phase="merge")
    with pytest.raises(Exception) as ei:  # fails later, on the stub file contents
        sg.run_ucscgb_chunked(a)
    assert "refusing to merge" not in str(ei.value), (
        f"a complete set of real chunk files was rejected as incomplete: {ei.value}"
    )


def test_merge_does_not_accept_details_left_uncompressed(tmp_path):
    """The merge seeks into BGZF and rebases byte offsets against it; a plain .tab is
    not a substitute, so the uncompressed name must not satisfy the requirement."""
    bed = _write_bed(tmp_path / "g.bed", 20)
    root = tmp_path / "chunks"
    root.mkdir()
    plain = tuple(f for f in REAL_CHUNK_FILES if not f.startswith("ascas12a_crisprDetails"))
    for i in range(2):
        _fake_chunk(str(root), i, files=plain + ("ascas12a_crisprDetails.tab",))

    a = _args(bed, tmp_path, ucscgb_phase="merge")
    with pytest.raises(RuntimeError, match="ascas12a_crisprDetails.tab.gz"):
        sg.run_ucscgb_chunked(a)


def test_merge_refuses_an_incomplete_set(tmp_path):
    """The dangerous failure of a sharded build: one shard never ran, the merge
    succeeds anyway, the bigBed is well-formed, and the percentiles are ranked over
    part of the genome. Nothing downstream would notice, so the merge has to."""
    bed = _write_bed(tmp_path / "g.bed", 50)  # 5 chunks
    root = tmp_path / "chunks"
    root.mkdir()
    for i in (0, 1, 2, 4):  # chunk 3 missing
        _fake_chunk(str(root), i)

    a = _args(bed, tmp_path, ucscgb_phase="merge")
    with pytest.raises(RuntimeError, match="refusing to merge 4 of 5 chunks"):
        sg.run_ucscgb_chunked(a)


def test_merge_refuses_a_chunk_missing_one_artifact(tmp_path):
    """A chunk dir that exists but lacks its bigBed is just as bad as an absent one."""
    bed = _write_bed(tmp_path / "g.bed", 30)  # 3 chunks
    root = tmp_path / "chunks"
    root.mkdir()
    for i in range(3):
        _fake_chunk(str(root), i)
    os.remove(os.path.join(sg._ucscgb_chunk_paths(str(root), 1)[1], "ascas12a.bb"))

    a = _args(bed, tmp_path, ucscgb_phase="merge")
    with pytest.raises(RuntimeError, match="missing \\['ascas12a.bb'\\]"):
        sg.run_ucscgb_chunked(a)


def test_merge_refuses_a_chunk_missing_its_guides_bed(tmp_path):
    """guides.bed carries the _offset column, the only link between a guide and its
    off-target details. A chunk without it cannot be merged at all."""
    bed = _write_bed(tmp_path / "g.bed", 30)
    root = tmp_path / "chunks"
    root.mkdir()
    for i in range(3):
        _fake_chunk(str(root), i)
    os.remove(os.path.join(sg._ucscgb_chunk_paths(str(root), 2)[1], "guides.bed"))

    a = _args(bed, tmp_path, ucscgb_phase="merge")
    with pytest.raises(RuntimeError, match="missing \\['guides.bed'\\]"):
        sg.run_ucscgb_chunked(a)


def test_merge_refuses_chunks_built_by_different_code(tmp_path):
    """All chunks must declare the same autoSql. Two different schemas means they were
    not all built by the same binary or snapshot, which is exactly how five published
    genomes ended up with pre-fix Cas12a scores."""
    bed = _write_bed(tmp_path / "g.bed", 30)
    root = tmp_path / "chunks"
    root.mkdir()
    _fake_chunk(str(root), 0)
    _fake_chunk(str(root), 1)
    _fake_chunk(str(root), 2, as_text='table T\n"t"\n( string chrom; "DIFFERENT" )\n')

    a = _args(bed, tmp_path, ucscgb_phase="merge")
    with pytest.raises(RuntimeError, match="chunks disagree about the schema"):
        sg.run_ucscgb_chunked(a)


def test_split_phase_stops_before_scoring(tmp_path):
    """--ucscgb_phase split writes the chunk BEDs and returns, scoring nothing."""
    bed = _write_bed(tmp_path / "g.bed", 25)
    a = _args(bed, tmp_path, ucscgb_phase="split")
    sg.run_ucscgb_chunked(a)
    root = tmp_path / "chunks"
    assert sorted(p.name for p in root.glob("chunk_*.bed")) == [
        "chunk_000.bed", "chunk_001.bed", "chunk_002.bed",
    ]
    assert not list(root.glob("track_*")), "split must not score a chunk"


def test_unknown_enzyme_is_rejected_before_any_merge_work(tmp_path):
    bed = _write_bed(tmp_path / "g.bed", 10)
    a = _args(bed, tmp_path, ucscgb_phase="merge", ucscgb_enzyme="nosuch")
    with pytest.raises(ValueError, match="unknown --ucscgb_enzyme"):
        sg.run_ucscgb_chunked(a)
