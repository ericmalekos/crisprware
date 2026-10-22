# SLURM drivers for a genome-wide `--ucscgb` track build

Everything here is meant to be run again for a new genome or a new enzyme. The one-off
half of the October 2026 AsCas12a rescan (watchers keyed to job IDs, retrofits the build
now does natively, the validation gate, the ISA diagnostics) is kept separately in the
`cas12a-rescan-campaign-2026-10` repository and should not be reached for.

## The scripts

| script | what it does |
|---|---|
| `make_snapshot.sh` | Freezes the live checkout as a snapshot a campaign can run from, and refuses to carry a binary forward unless `git diff <base>..HEAD -- crispr-ots/` is empty |
| `submit_shard_chain.sh` | Submits one genome as split then a shard array then a merge |
| `pipeline_shard.sh` | The per-phase worker the chain submits, `PHASE=split\|shard\|merge` |
| `verify_track.sh` | Verifies a finished track, including the row count no other check makes |
| `compare_tracks.py` | Compares two builds of one track on content, with columns graded hard and soft |
| `join_details.py` | Resolves `_offset` by merge join so two builds' details files can be compared |
| `push_ascas12a_tracks.sh` | rsyncs finished tracks to hgwdev. **Read the warning at the top of it before running.** |

## Normal use

```bash
# once per code change, from the live checkout
./make_snapshot.sh <base-snapshot-sha>

# once per genome
./submit_shard_chain.sh <org> <chunk_max> <shard_hours> gpu

# once the merge lands
./verify_track.sh <org> <suffix>
```

`PROJ` sets the campaign data directory and `REPO` the frozen code snapshot; both have
defaults and both are overridable, so a later campaign does not edit a script that a
running array may be reading. `submit_shard_chain.sh` resolves `pipeline_shard.sh` next
to itself rather than inside `PROJ`, so the pair travels together when the data moves.

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
