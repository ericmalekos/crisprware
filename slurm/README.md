# SLURM drivers for a genome-wide `--ucscgb` track build

Everything here is meant to be run again for a new genome or a new enzyme. The one-off
half of the October 2026 AsCas12a rescan (watchers keyed to job IDs, retrofits the build
now does natively, the validation gate, the ISA diagnostics) is kept separately in the
`cas12a-rescan-campaign-2026-10` repository and should not be reached for.

## The scripts

| script | what it does |
|---|---|
| `genomes.sh` | The per-genome mapping (db, tree, stem, chrom.sizes, chunk size, index and guide paths) in ONE place, sourced by the three scripts below |
| `make_snapshot.sh` | Freezes the live checkout as a snapshot a campaign can run from, and refuses to carry a binary forward unless `git diff <base>..HEAD -- crispr-ots/` is empty |
| `build_index.sh` | Builds the **TTTV** off-target index for one genome and verifies the PAM set it actually wrote |
| `submit_shard_chain.sh` | Submits one genome as split then a shard array then a merge |
| `pipeline_shard.sh` | The per-phase worker the chain submits, `PHASE=split\|shard\|merge` |
| `verify_track.sh` | Verifies a finished track, including the row count no other check makes |
| `compare_tracks.py` | Compares two builds of one track on content, with columns graded hard and soft |
| `join_details.py` | Resolves `_offset` by merge join so two builds' details files can be compared |
| `push_ascas12a_tracks.sh` | rsyncs finished tracks to hgwdev. **Read the warning at the top of it before running**: it names two earlier track generations that must never be pushed. |

## Normal use

```bash
# once per code change, from the live checkout
./make_snapshot.sh <base-snapshot-sha>

# once per genome, the TTTV off-target index (CPU, seconds to minutes, no GPU)
sbatch -J idx_<db> -p gpu -c 4 --mem=48G -t 01:00:00 --export=ALL,ORG=<org> build_index.sh

# once per genome, the scan
./submit_shard_chain.sh <org> <chunk_max> <shard_hours> gpu

# once the merge lands
./verify_track.sh <org> <suffix>
```

`PROJ` sets the campaign data directory, `REPO` the frozen code snapshot and `SUFFIX`
the track directory name (`track_<db>_cas12a_$SUFFIX`, default `tttv`); all have
defaults and all are overridable, so a later campaign does not edit a script that a
running array may be reading. `submit_shard_chain.sh` resolves `pipeline_shard.sh` and
`genomes.sh` next to itself rather than inside `PROJ`, so they travel together when the
data moves.

Useful env vars on the chain submitter: `SHARDS_ONLY=1` to skip the split for a genome
already split, `RESUME=1` to submit only the chunks that are not already complete,
`THROTTLE=<n>` to cap concurrent array tasks, and `MERGE_PART` / `MERGE_HOURS` for the
merge's placement.

## What this shape is for

The chunk work was always independent: each chunk is its own process with its own output
directory and no shared state. But the in-process worker count is `max(1, len(gpus))`,
which is 1 on the CPU scanner, so a 68-chunk genome ran strictly serially. The array
runs the same chunks concurrently.

Two properties matter before using it. The split is the only writer of `chunk_*.bed`, so
it must finish before any shard starts; N shards splitting concurrently would interleave
their writes and produce chunks that are individually well-formed and collectively
wrong. And the campaign is resumable at shard granularity, because chunks are kept: a
shard that fails costs that shard rather than the genome, and the merge refuses an
incomplete set anyway.

## Findings that shaped these scripts

Each of these cost real time to learn and is encoded in the code rather than left to
whoever runs it next.

**Request `--gres=gpu:1`, never a model.** Naming a model restricts the job to a subset
of nodes, and it then queues behind that subset while other GPUs sit idle. On 2026-10-17
a 4,760-task array pinned to `A5500` got 3 concurrent slots instead of 24, turning about
25 task-hours into a 16 hour wall time while 8 A100s sat unreachable to it.

**The merge belongs wherever there are free CPUs, which is usually not `medium`.** On
prism, `short`, `medium` and `long` share one CPU pool and it is normally saturated,
while the GPU nodes carry a separate pool with hundreds of idle CPUs. The merge is pure
I/O: it reads each chunk's `guides.bed` and BGZF details, rebases offsets, re-ranks
percentiles and writes the bigBed. No GPU, no TensorFlow, and therefore no avx2/avx512
sensitivity, which is what makes any partition safe for it. Moving one merge from
`medium` to `gpu` with no `--gres` took it from a 17 hour queue estimate to starting in
40 seconds.

**Size time limits from measurement, not from caution.** An oversized request has to
find an oversized window. Measured on the eight-genome campaign: shards 7 to 27 minutes
at `chunk_max` 1,000,000, merges 23 to 53 minutes for 37 to 142 chunks. The 24 hour
per-shard default this work started with was two orders of magnitude off.

**`PriorityFlags` is empty on this cluster**, so the age factor accrues from eligible
time rather than submit time. A job blocked on a dependency or held banks zero age,
which is why `SHARDS_ONLY=1` exists: an array submitted with
`--dependency=afterok:<split>` starts its priority clock only when the split finishes,
and one genome waited 18.3 hours on `Reason=Priority` before placing.

**Pin the TensorFlow thread count.** Cleavage scores depend both on the host ISA group
and on intra-op thread count. avx2+fma GPU nodes agree byte for byte; avx512 nodes
differ on about 53% of values by up to 2.3e-05. `pipeline_shard.sh` sets
`TF_NUM_INTRAOP_THREADS` from `SLURM_CPUS_PER_TASK` rather than from `nproc`, which does
not respect the cgroup on every node.

**Assert which code is running.** `python -m` puts the current directory at `sys.path[0]`
ahead of `PYTHONPATH`, so a job can silently import a live file instead of the frozen
snapshot. `pipeline_shard.sh` resolves `crisprware.score_guides.__file__` and fails if it
is not inside the snapshot. It also checks the scanner binary's md5 against an
allow/deny list, because the campaign exists only because a pre-fix binary was carried
forward through five successive snapshots while the Python was correctly re-frozen each
time.

**Verify row counts, not chunk counts.** The merge checks that every expected chunk is
present, which cannot detect rows dropped inside a chunk. `verify_track.sh` compares the
bigBed's `itemCount` against the source guide BED, which can.

**Index AsCas12a on TTTV, and make the code check it rather than assume it.** AsCas12a
does not cleave a TTTT-prefixed protospacer, so a TTTT site in the off-target index is a
false off-target. The scanner excludes TTTT from exactly one of its outputs,
`specificity_tttv`; the mismatch buckets behind the close-off-target grey, the off-target
list in the details popup, the off-target count, max CFD, the saturation cap and the 0-mm
screen behind `unique` all count every site they are handed. Eight genomes were indexed
on TTTN while three separate comments in `ucsc_track.py` asserted "TTTV-only index", and
nothing could tell: on hg38 chr1 that greyed 12.49% of guides as hard to target whose own
published specificity said otherwise, and TTTT turned out to be 37.5% of a TTTN index
rather than the naive quarter, because poly-T runs concentrate in soft-masked repeats.
Two things follow. Fix it at the index, not at each consumer, because eight separate
guards are eight chances to forget one. And read the assumption off the file:
`crisprware/utils/crot_index.py` parses the `.crot` header and
`score_guides --ucscgb_enzyme ascas12a` refuses a TTTT-bearing index before scanning,
which is a check `build_index.sh` repeats on what it just wrote.

### Removed binary: `crispr-ots-gpu` md5 `22e7030f06b1`

Deleted 2026-10-23 from `crispr-ots/target/x86_64-unknown-linux-gnu/release/`. It was a
2026-06-06 build, predating `ad3ecc6` ("crispr-score: fix Cas12a off-target CFD matrix
application (complement the d-base)"), so its Cas12a off-target CFDs are wrong. It stays
on the quarantine md5 list in `pipeline_shard.sh` and `build_index.sh` so a restored copy
is still refused; the list is the guard, not the file's absence.

Full md5 `22e7030f06b1e994a6c81f0ab3d21807`, 9,744,304 bytes, mtime 2026-06-06 16:33:55.
It was harmless for SpCas9 (`ad3ecc6` touches `cas12a.rs`, not `cfd.rs`, and is the only
`crispr-score` commit between that build and 2026-08-09), but nothing should be able to
reach it for Cas12a. Two references outside this repo pointed at it:
`crispr-ots-multigenome/config/params.env` (`OTS_GPU`) and
`crispr-ots-multigenome/bin/make_tutorial.py:145`, which emitted a tutorial command line
invoking it.
