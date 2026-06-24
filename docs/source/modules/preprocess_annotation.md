# preprocess_annotation

Collapse a multi-isoform GTF/GFF into simple per-gene **models**, optionally filtering isoforms by
RNA-seq expression, so the later steps target one well-defined region per gene.

```bash
crisprware preprocess_annotation -g tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf \
  -m metagene consensus longest shortest
```

```text
chrIII_ce11.ncbiRefSeq/
  chrIII_ce11.ncbiRefSeq_meta.gtf            # metagene: union of all exons
  chrIII_ce11.ncbiRefSeq_consensus.gtf       # exons shared by every isoform
  chrIII_ce11.ncbiRefSeq_longestCDS.gtf      # longest-CDS isoform per gene
  chrIII_ce11.ncbiRefSeq_shortestCDS.gtf     # shortest-CDS isoform per gene
  chrIII_ce11.ncbiRefSeq_genes_without_consensus_model.txt
```

## Options you reach for

- `-m/--model` `metagene consensus longest shortest` -- which models to emit (any subset).
- **RNA-seq filter**: `-t/--tpm_files quant.sf ...` (Salmon / Kallisto / FLAIR), then `--median`,
  `--mean`, `--min`, `--max` cutoffs, or `-n/--top_n` to keep the N highest-expressed isoforms per gene.
- **dCas targeting**: `--tss_window <up> <down>` / `--tes_window <up> <down>` emit BED windows around
  the TSS / TES for CRISPRi/a.

```{note}
Feed the model you choose to `generate_guides -k` (restrict guides to it) and to `rank_guides -t`
(group output per gene/transcript). Make sure the GTF and any TPM files share transcript IDs.
```
