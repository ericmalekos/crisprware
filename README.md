# CRISPRware

[![Tests](https://github.com/ericmalekos/crisprware/actions/workflows/tests.yml/badge.svg)](https://github.com/ericmalekos/crisprware/actions/workflows/tests.yml)
[![Lint](https://github.com/ericmalekos/crisprware/actions/workflows/lint.yml/badge.svg)](https://github.com/ericmalekos/crisprware/actions/workflows/lint.yml)

CRISPRware is a toolkit to preprocess NGS data and identify, score, and rank guide RNAs (gRNAs) for
CRISPR experiments. It supports RNA-seq, Ribo-seq, ATAC-seq, DNase-seq, ChIP-seq, and other genomic
preprocessing.

**v0.2**: 12 in-tree on-target scorers spanning SpCas9 and Cas12a (all bit-equivalence-tested against
their upstream references), plus an in-tree Rust off-target backend (`crispr-ots`, ~6x faster than the
previous `guidescan` shell-out) with bundled Cas12a CFD matrices. `--cas9_scorer` and `--cas12a_scorer`
accept multiple values, so one call can emit any combination of scores.

## Documentation

Full documentation lives at **[crisprware.readthedocs.io](https://crisprware.readthedocs.io)**:
[installation](https://crisprware.readthedocs.io/en/latest/installation.html), a page per module,
[scoring methods](https://crisprware.readthedocs.io/en/latest/scoring/index.html),
[NGS-guided design](https://crisprware.readthedocs.io/en/latest/ngs.html), the
[Cas12a UCSC track pipeline](https://crisprware.readthedocs.io/en/latest/ucsc_cas12a/index.html), and
the full [command reference](https://crisprware.readthedocs.io/en/latest/command_reference.html).

## Implemented scoring methods

On- and off-target methods built into `score_guides`; each links to its
[documentation](https://crisprware.readthedocs.io) page:

| Method | Type | Nuclease | Reference |
|--------|------|----------|-----------|
| [RuleSet3](https://crisprware.readthedocs.io/en/latest/scoring/ontarget.html#ruleset3) | on-target | SpCas9 | DeWeirdt 2022 |
| [DeepSpCas9](https://crisprware.readthedocs.io/en/latest/scoring/ontarget.html#deepspcas9) | on-target | SpCas9 | Kim 2019 |
| [DeepHF](https://crisprware.readthedocs.io/en/latest/scoring/ontarget.html#deephf) (3 variants) | on-target | SpCas9 | Wang 2019 |
| [DeepCpf1](https://crisprware.readthedocs.io/en/latest/scoring/ontarget.html#deepcpf1) | on-target | Cas12a | Kim 2018 |
| [enPAM+GB](https://crisprware.readthedocs.io/en/latest/scoring/ontarget.html#enpam-gb) | on-target | enAsCas12a | DeWeirdt 2021 |
| [enseq-DeepCpf1](https://crisprware.readthedocs.io/en/latest/scoring/ontarget.html#enseq-deepcpf1) | on-target | Cas12a | Chen 2025 |
| [seq-DeepCpf1variants](https://crisprware.readthedocs.io/en/latest/scoring/ontarget.html#seq-deepcpf1variants) | on-target | 23 Cas12a variants | Chen 2025 |
| [crispr-ots](https://crisprware.readthedocs.io/en/latest/scoring/offtarget.html#crispr-ots) (CFD) | off-target | SpCas9 + Cas12a | DeWeirdt 2021 / Doench 2016 |
| [Guidescan2](https://crisprware.readthedocs.io/en/latest/scoring/offtarget.html#guidescan2) | off-target | SpCas9 + Cas12a | Schmidt 2023 |

## Quickstart

Install per the [installation docs](https://crisprware.readthedocs.io/en/latest/installation.html)
(Python env + the `crispr-ots` binary on your `PATH`).

### Inputs
- **FASTA** (required).
- **BED** (optional): regions of interest, to limit the gRNA search space.
- **GTF/GFF** (optional): gene annotation, to filter gRNAs by feature (e.g. exon, CDS).

### Walkthrough

Demonstrated on the ce11 chromosome III FASTA + NCBI GTF in `tests/test_data/ce11/` (the example
off-target index is limited to chrIII, not the full ce11 genome).

`index_genome` builds a **crispr-ots** off-target index (fast, exact-mismatch) into a
`<name>_crisprots/` directory. `score_guides` auto-detects the engine from each `-i` index, so it can
also score a Guidescan2 index (built with Guidescan2's own tooling) for RNA/DNA-bulge searches.

```
crisprware index_genome -f tests/test_data/ce11/chrIII_sequence.fasta -p NGG -l 20
```

Build gene models from the NCBI GTF:

```
crisprware preprocess_annotation -g tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf \
-m metagene consensus longest shortest
```

Default settings generate NGG protospacer guides:

```
crisprware generate_guides -f tests/test_data/ce11/chrIII_sequence.fasta \
-k tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf \
--feature CDS
```

Scoring takes ~5 minutes and uses 8 threads by default (`--threads <int>`). `--tracr` is `Chen2013`,
`Hsu2013`, or `both` (see [RuleSet3](https://github.com/gpp-rnd/rs3)). Pass more than one `-i` index to
get a separate `specificity_<index>` column from each engine in one run.

```
crisprware score_guides -b chrIII_sequence_gRNA/chrIII_sequence_gRNA.bed \
-i chrIII_sequence_crisprots/chrIII_sequence_crisprots --tracr Chen2013 --threads 8
```

Rank on scoring columns: `-c` paired with `-m` filters (here `RS3_score_Chen2013 >= 0` and
`specificity_chrIII_sequence_crisprots >= 0.2`); `-p 5 65 -f CDS` keeps the 5th-65th CDS percentile;
`--output_all` writes a TSV + histograms per filter stage.

```
crisprware rank_guides \
-k chrIII_sequence_scoredgRNA/chrIII_sequence_scoredgRNA.bed \
-t tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf \
-f CDS \
-c RS3_score_Chen2013 specificity_chrIII_sequence_crisprots \
-m 0 0.2 \
-p 5 65 \
-r RS3_score_Chen2013 \
--output_all
```
