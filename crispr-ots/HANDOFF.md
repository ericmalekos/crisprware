# crispr-ots — agent handoff

You are picking up `crispr-ots`, a Rust off-target scanner that lives
inside the `crisprware` Python project. This document is the briefing —
read it end-to-end before touching code. There is **no Claude memory
store on this machine**; everything you need to know about decisions
made on the dev branch is in here or in the linked file paths.

---

## What this project is

`crisprware` is a Python pipeline for designing CRISPR sgRNAs. Its
core stages:

```
preprocess_annotation → generate_guides → index_genome → score_guides → rank_guides
                                              │              │
                                              │              ├── RS3 cleavage score (lightgbm)
                                              │              └── off-target specificity score  ← THIS PROJECT
                                              │
                                              └── builds an off-target index ─────────────────┘
```

Historically `index_genome` and `score_guides` shelled out to
**GuideScan2** (`guidescan index` / `guidescan enumerate`) for the
off-target side. We replaced both with `crispr-ots`, a Rust workspace
that lives at `crispr-ots/` inside this repo. It is **faster than
GuideScan2 and FlashFry on every dimension** (~13× faster enumerate
on mouse, ~5× faster index build) and produces byte-compatible output
so the rest of crisprware's pipeline doesn't change.

The branch you're on is `dev/crispr-ots`. It is **N commits ahead of
`main`**; the last commit when this doc was written is `ebfaa63`
("environment.yml: swap guidescan for the conda-built crispr-ots").
Plan is to keep iterating here and merge to `main` once the cluster
benchmarks back up the perf claims and a couple of follow-ups land
(see "Open work" below).

---

## Repository layout

```
crisprware/                         the Python package (unchanged from upstream
  crisprware/                       except for two subprocess.run() lines —
    index_genome.py:85              calls `crispr-ots index` instead of `guidescan index`
    score_guides.py:186             calls `crispr-ots enumerate` instead of `guidescan enumerate`
  environment.yml                   conda env spec; depends on `crispr-ots` from a local channel
  conda-recipe/                     conda-build recipe for crispr-ots
    crispr-ots/meta.yaml            the recipe
    build_local_channel.sh          driver: builds `./conda-channel/`
  parasol_scripts/                  legacy / utility Python scripts
    score_flashfry_cfd.py           the Cas12a CFD reference impl we ported to Rust
    off_targ_2xNLS_Cas12a.csv       Cas12a activity matrix (23 positions × 12 mismatch types)
    off_targ_enCas12a.csv           alt Cas12a matrix (engineered broad-PAM variant)
  crispr-ots/                       the Rust workspace (this directory)
    Cargo.toml                      5-crate workspace
    .cargo/config.toml              sets target-cpu=x86-64-v3 (Haswell+)
    BENCHMARKS.md                   ★ READ THIS; every perf claim in this doc traces back to it
    NOTES.md                        ★ READ THIS; caveats from FlashFry/GuideScan2/Cas12a science
    HANDOFF.md                      ← this file
    crates/
      crispr-encoding/              2-bit DNA, Site (2 × u64, up to 28 bp), Enzyme presets, mismatch primitive
      crispr-db/                    BinTable (build) + MmapDb (query, zero-copy), v2 format, BinSource trait
      crispr-scan/                  BinScanner: bin-prefix prefilter + SIMD popcount inner loop (rayon + wide)
      crispr-score/                 SpCas9 CFD (Cfd) + Cas12a CFD-like (Cas12aCfd); matrices in data/
      crispr-cli/                   CLI: `index` (alias `build`) + `enumerate` (alias `discover`)
    mm39-bench/                     mouse-genome benchmark drivers + 1000-guide fixtures
      make_guides.py                generates 1000 random SpCas9 NGG guides (seed 42)
      make_guides_cpf1.py           generates 1000 random Cas12a TTTN guides (seed 42)
      binwidth_sweep.sh             bin-width sweep driver
      run_flashfry_chunked.sh       FlashFry 8-process discover wrapper
      compare_specs.py              cross-tool specificity diff
      random_1000.fasta             SpCas9 fixture (committed)
      random_1000.kmers.csv         SpCas9 GS2-format fixture (committed)
      random_1000_cpf1.fasta        Cas12a fixture (you'll need to regenerate after pulling — see below)
      random_1000_cpf1.kmers.csv    Cas12a GS2-format fixture (same)
```

The mouse `.crot` indexes (4-6 GB each) are **not committed** and need
to be rebuilt on this machine. See "How to rebuild fixtures" below.

---

## What's in the engine right now (feature surface)

```
crispr-ots index FASTA --index PREFIX [--enzyme NAME] [--bin-width 1..15]
```

- Builds an mmap-friendly v2-format off-target database. Emits
  `<PREFIX>.crot` plus zero-byte stubs at `<PREFIX>.{forward,reverse,gs}`
  so crisprware's `check_files_exist` (which looks for those three
  GuideScan2 names) succeeds unmodified.
- Enzymes: `spcas9ngg` (default, 20-nt protospacer + 3-nt NGG PAM),
  `spcas9nag` (same but NAG), `spcas9` (NGG ∪ NAG), `cpf1` /
  `cas12a` (4-nt TTTN PAM + **23-nt protospacer** — see "Cas12a
  protospacer" below).
- Bin-width is the inner-loop prefilter sharpness; default 7 (good for
  chr-scale), empirical sweet spot is **11 for full mouse / human**.

```
crispr-ots enumerate INDEX_PREFIX
    --kmers-file CSV | --queries FASTA
    [--mismatches 4] [--threads 0]
    [--format csv|tsv|both]
    [--spec-convention guidescan|flashfry]
    [--score cfd | cfd-cas12a | cfd-cas12a:encas12a]
    [--threshold N | -1] [--max-off-targets-per-bin N | -1]
    [--rna-bulges 0] [--dna-bulges 0]    # accepted but reject non-zero
    [--alt-pam ... ] [--mode ...] [--max-off-targets ...]    # accepted, GS2-compat noop
    --output PATH
```

- Input: crisprware/GuideScan2 CSV (preferred for crisprware
  integration) or a FASTA. CSV format is documented in `kmers_csv.rs`.
- Output formats:
    - `csv` — one row per (guide, hit). The GuideScan2 schema
      `id,sequence,match_chrm,match_position,match_strand,
      match_distance,specificity`. This is what `score_guides.py`
      reads back.
    - `tsv` — one row per guide; FlashFry-style with `otCount`,
      `otSequences`, and optional CFD columns. Useful for the
      "report all off-targets" use case.
    - `both` — CSV at `--output`, TSV at `<output>.detail.tsv`.
- Specificity convention selectable:
    - `--spec-convention flashfry` → `1 / (1 + Σ_OT cfd × count)`.
      Default for `tsv`.
    - `--spec-convention guidescan` → `1 / Σ_all_rows cfd`. Default
      for `csv` (matches what crisprware expects).
- Scoring:
    - `--score cfd` — Doench 2016 CFD for SpCas9 only.
    - `--score cfd-cas12a` (alias `cfd-cas12a:2xnls`) — Cas12a
      CFD-like with the `2xNLS-Cas12a` activity matrix.
    - `--score cfd-cas12a:encas12a` — same model, engineered-variant
      matrix. Both ship as embedded CSVs in `crispr-score/data/`.
    - Emits `cas12a_max`, `cas12a_spec_tttn`, `cas12a_spec_tttv`
      columns. TTTV variant excludes off-targets whose PAM is TTTT
      (Cas12a tolerates those poorly).
- `--threshold N` — guidescan-compat filter; drop guides whose nearest
  *non-on-target* hit is within N mismatches. **Includes mm=0
  multi-mappers** (a guide whose sequence appears at ≥ 2 genomic
  positions at mm=0 is dropped at any non-negative threshold).
  Behavior fixed by commit `c3c7280` after a 7-guide diff vs real
  GuideScan2 on the ce11/chrIII smoke test.
- `--max-off-targets-per-bin N` (default **500**) — cap distinct
  off-target sequences per `(guide, mismatch-count)` bin. Critical
  for low-complexity guides which would otherwise emit 800 k+
  off-targets at mm=4 on mouse. Set `-1` to disable.

Read `crates/crispr-cli/src/main.rs` and `discover.rs` for the
authoritative arg parsing.

---

## What is *not* implemented yet (intentional)

- **DNA/RNA bulges**. `--rna-bulges N` / `--dna-bulges N` accept N=0
  only; non-zero rejected at runtime. GuideScan2 has bulge support
  via FM-index recursion; we don't (the bin-scan engine is
  mismatch-only by design). To add: would need a second `Scanner`
  trait impl behind an FM-index. Out of scope for now.
- **GS2 `--mode succinct|complete`**. Accepted as a no-op; we always
  enumerate every off-target up to `--mismatches` (i.e. always
  "complete").
- **GS2 `--alt-pam`**. Accepted as a no-op. We match the enzyme's
  IUPAC-expanded PAM set during the scan; passing `--alt-pam NGG`
  isn't necessary the way it is for GuideScan2.
- **GPU scan path**. Analyzed but not built — see "GPU feasibility
  analysis" in BENCHMARKS.md. Short version: the 4 GiB GTX 1650 Ti on
  the original dev box is too small for the mouse `.crot`; an A100 /
  A5500 has the headroom but the prototype hasn't been written. **This
  is the highest-impact open work item on the new hardware.**

---

## Correctness pinning (your safety net)

`cargo test --workspace --release` runs ~143 tests in 5-10 seconds.
The ones to watch carefully when you change anything:

| Test | What it locks down |
|---|---|
| `crispr-score::cas12a_parasol_validate` | Per-pair Cas12a CFD matches the parasol-script reference to ≤ 4e-11 across 5 hand-picked cases including a PAM-distal 4-mismatch case at matrix positions 21..23. **Do not break this.** |
| `crispr-cli::cli_cfd_chr22_1000::cfd_and_otcount_match_flashfry_on_1000_random_chr22_guides` | SpCas9 CFD bit-equivalence with FlashFry across 930 chr22 guides, max \|Δ spec\| = 1 ULP. Skips if `/tmp/ffry-quickstart/` is missing — regenerate fixtures if you want this test to run. |
| `crispr-cli::cli_guidescan2_equivalence` | Off-target position-count + spec match GuideScan2 to its 6-sigfig print precision floor. |
| `discover::tests::threshold_drops_multimapping_on_target` | The mm=0-multi-mapper threshold filter (the gotcha from `c3c7280`). |
| `discover::tests::max_per_bin_caps_distinct_off_target_sequences` | The cap behaves correctly per `(guide, mm)` bin. |
| `cas12a::tests::*` | Cas12a scorer parses both bundled matrices, multi-mismatch product is right, TTTV excludes TTTT. |

If you change scoring, the bin layout, or the threshold/cap logic and
*any* of those go red, stop and figure out why.

`cargo clippy --workspace --release --tests -- -D warnings` is also
green; keep it that way.

---

## Open work items (in priority order)

### 1. Run the mouse benchmark on the cluster's actual CPUs

Numbers in `BENCHMARKS.md` are from an i7-10870H laptop (Comet Lake,
AVX2 but no AVX-512 VPOPCNTDQ). The cluster's host CPUs (paired with
the A100 / A5500) are almost certainly newer (Ice Lake-SP / Sapphire
Rapids / Zen 4) and have **VPOPCNTDQ**. Even without a code change,
just running on hardware with more cores + faster memory should shift
the headline numbers.

**Concrete:** copy the mm39 fixture or rebuild it, run
`crispr-ots index --enzyme spcas9ngg --bin-width 11 --index mm39
GRCm39.fa.gz`, then `crispr-ots enumerate --threads $(nproc) ...` and
append the rows to `BENCHMARKS.md`.

### 2. AVX-512 VPOPCNTDQ scan path

Today's SIMD scan uses `wide::u64x4` + SWAR popcount because the dev
box doesn't have hardware lane-wise popcount. On a VPOPCNTDQ-capable
CPU the SWAR can be replaced by a single intrinsic per lane.

**Concrete:** add `multiversion` to `crispr-scan/Cargo.toml`, write a
`#[multiversion(targets = "simd")]` variant of `pop_pairs_x4` that
uses `_mm512_popcnt_epi64` on Ice Lake+ and falls back to the current
SWAR otherwise. The hot function is in
`crates/crispr-scan/src/bin_scanner.rs:pop_pairs_x4`. Expected gain:
3-5× on the scan kernel for VPOPCNTDQ-capable cores, ~1.5-2× overall.

### 3. GPU prototype on the A5500 / A100

Now that the hardware constraint is gone (24 GB on A5500, 40 or 80 GB
on A100, all of which comfortably hold a mouse `.crot`), the GPU path
becomes worth building. The feasibility analysis in BENCHMARKS.md
estimates 3-5× overall speedup on mouse-scale workloads.

**Concrete:** prototype in `cudarc` (NVIDIA-only, lowest overhead)
or `wgpu` (cross-platform, slightly more overhead). The kernel is
the mismatch-popcount inner loop. Suggested approach:

1. Build a stripped-down crate `crispr-scan-gpu` exposing
   `Scanner` trait impl that copies bin entries + guides to GPU,
   runs a kernel, copies hits back.
2. Start on chr22 (small fixture, fits in any VRAM) for a quick CPU vs
   GPU comparison. If it wins there, scale to mouse.
3. Keep the on-disk `.crot` format unchanged; the GPU just reads bin
   chunks streamed from the mmap'd CPU side.
4. Pin perf gains in a new `BENCHMARKS.md` section.

If the prototype shows ≥ 3× wall-time gain on mouse, it's worth
landing. If it shows < 2×, document the result and shelve. Don't
fight the hardware.

### 4. Bioconda submission

`conda-recipe/crispr-ots/meta.yaml` is already written and validated
locally. Once the API stabilizes (after either of the two perf items
above lands or doesn't), open a bioconda-recipes PR. Then
`environment.yml` drops the `${CRISPR_OTS_CHANNEL}` hack and just
lists `- crispr-ots` under `bioconda::`.

### 5. Smaller things

- **Zstd-compressed `.crot` variant.** Cold storage of mouse indexes
  is 4-6 GB; zstd should cut that to ≤ 1 GB at the cost of a
  decompress on open. Useful for distribution. Behind a feature
  flag so the mmap-friendly variant stays the default.
- **Per-bin cap inside the scanner**. The current `--max-off-targets-per-bin`
  trims at the post-scan grouping stage; the scan still pays for
  finding the dropped entries. Pushing the cap into the bin scanner
  (per-thread `(guide_idx, mm)` counters, skip-emit when at cap)
  would give real scan-time savings on low-complexity guides, not
  just bounded output. ~5-15× speedup possible on overflow guides.
- **23-nt-protospacer SpCas9 enzyme.** Cas12a is now 23-nt; SpCas9 is
  still 20. The encoding supports any length up to 28. Not biologically
  motivated (the SpCas9 CFD matrix is 20 positions) but the engine is
  uniform.

---

## Non-obvious decisions, in case you're tempted to change them

- **23-nt Cas12a protospacer**. We diverged from FlashFry's 20-nt
  Cpf1ParameterPack here. The 23-nt enzyme matches the published
  Cas12a activity matrices (which have 23 positions, not 20). The
  2× scan speedup we measured after the switch is a bonus from a
  sharper bin prefilter. Don't roll this back.
- **Threshold filter includes mm=0 multi-mappers**. If a guide's
  exact protospacer appears at 2+ genomic positions, ALL copies of
  that guide get filtered. This matches guidescan2's `--threshold`
  semantics; we discovered it via the ce11/chrIII smoke test (7-guide
  diff). Lives in `discover.rs:225-247` with a long comment.
- **String-equality match logic in Cas12a CFD**. The bundled matrices
  contain entries like `rA:dA` that LOOK like Watson-Crick same-base
  mismatches but are dead code under the upstream parasol script's
  `tb != ob → lookup, else skip` convention. We follow that convention
  verbatim so output matches the parasol pipeline byte-for-byte.
  Don't try to "fix" this without changing the upstream too.
- **target-cpu=x86-64-v3 in `.cargo/config.toml`**. Haswell+ baseline.
  Required for the AVX2 codegen in `wide::u64x4`. If you need to
  support older CPUs, replace with `multiversion` and lose the
  config-level constraint. **All modern HPC compute nodes satisfy v3.**
- **GS2-compat args are accepted-but-noop**. `--alt-pam`, `--mode`,
  `--max-off-targets`, `--rna-bulges 0`, `--dna-bulges 0`. We accept
  them so crisprware's existing `subprocess.run` invocation works
  unchanged; they don't do anything inside the engine. **Don't
  delete these — `crisprware/score_guides.py:186` passes them.**

---

## How to rebuild data on this machine

If `mm39-bench/` is missing indexes or the SpCas9 random_1000 fixture
is stale, here's the order:

```bash
# 1. Reference genome (~3 GB decompressed)
gunzip -k -c GRCm39.primary_assembly.genome.fa.gz > mm39-bench/GRCm39.fa

# 2. SpCas9 fixture (already committed but if you need to regenerate)
cd mm39-bench
python3 make_guides.py GRCm39.fa                  # 1000 guides, seed 42, ~5 min

# 3. Cas12a fixture (need this if it's not committed or if it's a
#    24-bp variant from before the 23-nt switch)
python3 make_guides_cpf1.py GRCm39.fa             # 1000 guides, seed 42

# 4. crispr-ots SpCas9 index (~5 min, 6.3 GB on disk)
crispr-ots index --enzyme spcas9ngg --bin-width 11 \
    --index mm39.crispr-ots GRCm39.primary_assembly.genome.fa.gz

# 5. crispr-ots Cas12a index (~5 min, 4.0 GB on disk)
crispr-ots index --enzyme cpf1 --bin-width 11 \
    --index mm39.crispr-ots.cpf1 GRCm39.primary_assembly.genome.fa.gz
```

For the FlashFry / GuideScan2 indexes referenced in BENCHMARKS.md,
see commands in BENCHMARKS.md `## Reproducing`. **You don't need
those for development; they're only for running cross-tool perf
comparisons.** The local cross-tool comparison was already done; only
re-run them if you're updating the benchmark numbers.

---

## How to ship a change

1. Edit code. `cargo test --workspace --release` must stay green.
2. `cargo clippy --workspace --release --tests -- -D warnings` must
   stay green.
3. If touching CLI surface or output schema: also run a smoke through
   the crisprware Python wrappers using
   `PYTHONPATH=$REPO_ROOT python -m crisprware.cli score_guides ...`
   on the `ce11/chrIII` fixture in `crisprware/tests/test_data/`.
4. If touching perf-critical code: re-run the relevant
   `mm39-bench/*.sh` driver, append a row to BENCHMARKS.md's snapshot
   history table.
5. `git commit` — descriptive title (≤ 70 chars), body explains
   *why* not *what*. **Do not add `Co-Authored-By: Claude` trailers**
   (the user has a global instruction against these on all repos).
6. **Do not push or merge without an explicit "push" instruction**
   from the user. `dev/crispr-ots` is the working branch; main is
   protected.

The first agent on a fresh clone should run the smoke test in step 3
once before changing anything, just to confirm the install is good.

---

## Final pointer

`BENCHMARKS.md` is the source of truth for every perf claim. `NOTES.md`
is the source of truth for every "why is it like this" gotcha. Read
them both before touching anything controversial. If you find
yourself wanting to make a 100-line change without consulting either,
stop and re-read this section.
