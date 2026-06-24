# Off-target scoring

Off-target specificity comes from scanning each guide against an index built by `index_genome`. The
engine is auto-detected per `-i` index, so passing one of each gives a `specificity_<index>` column
from each (0-1, where 1.0 = no off-targets).

## crispr-ots

Default engine (`index_genome --indexer crispr-ots`). Rust scanner, **exact mismatches only** (no
bulges), ~6x faster than Guidescan2. Emits a CFD-style aggregate specificity, with bundled matrices
for **SpCas9** (CFD) and **Cas12a** (`enCas12a` and `2xNLS-Cas12a`). The Cas12a CFD matrices are from
[DeWeirdt et al. 2021](https://www.nature.com/articles/s41587-020-0600-6)
([code](https://github.com/PeterDeWeirdt/cas12a_manuscript)); SpCas9 CFD from
[Doench et al. 2016](https://www.nature.com/articles/nbt.3437).

## Guidescan2

`index_genome --indexer guidescan2`. Trie-based; supports **RNA/DNA bulges**
(`score_guides --rna_bulges` / `--dna_bulges`) in addition to mismatches, at higher runtime/memory.
[Schmidt et al. 2023](https://www.biorxiv.org/content/10.1101/2022.05.02.490368v1). Prebuilt indices
for some genomes are available at <https://guidescan.com/downloads>.

## Controls (in `score_guides`)

- `--mismatches N` -- search radius.
- `--threshold N` -- discard a guide if any off-target sits within N mismatches (`-1` keeps all).
- `--rna_bulges` / `--dna_bulges` -- guidescan2 indices only (crispr-ots requires 0).
- `--alt_pams NAG ...` -- additional PAMs to consider as off-targets.
