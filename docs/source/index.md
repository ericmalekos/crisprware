# CRISPRware

CRISPRware identifies, scores, and ranks CRISPR guide RNAs (gRNAs) from a genome, with optional
NGS-guided targeting (RNA-seq, Ribo-seq, ATAC-seq, ChIP-seq). On-target activity is predicted by a
suite of deep-learning models for **SpCas9** and **Cas12a**; off-target specificity is scored by the
in-tree Rust scanner **`crispr-ots`** (or Guidescan2).

## Pipeline

```{mermaid}
flowchart LR
    A["<b>preprocess_annotation</b><br/><small>gene models, RNA-seq filter</small>"]
    I["<b>index_genome</b><br/><small>off-target index</small>"]
    G["<b>generate_guides</b><br/><small>find protospacers</small>"]
    S["<b>score_guides</b><br/><small>on-target + off-target</small>"]
    R["<b>rank_guides</b><br/><small>filter + select per target</small>"]
    A --> G --> S --> R
    I --> S
    classDef s fill:#eaf2fb,stroke:#3a6ea5,color:#16314a;
    class A,I,G,S,R s;
```

Each step is a `crisprware <module>` subcommand; chain them or run any in isolation. The pages below
walk each module with a runnable example on the bundled *C. elegans* chrIII test data.

```{toctree}
:caption: Getting started
:maxdepth: 1

installation
quickstart
tutorials
```

```{toctree}
:caption: Modules
:maxdepth: 1

modules/index
ngs
```

```{toctree}
:caption: Scoring
:maxdepth: 1

scoring/index
```

```{toctree}
:caption: Cas12a UCSC tracks
:maxdepth: 2

ucsc_cas12a/index
reading_the_track/index
```

```{toctree}
:caption: Reference
:maxdepth: 1

command_reference
references
```
