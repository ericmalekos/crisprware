# NGS-guided design

Restrict or prioritize guides using your own NGS data: expressed isoforms (RNA-seq), translated ORFs
(Ribo-seq), or accessible/bound regions (ATAC/ChIP, via BED windows).

## RNA-seq: keep expressed isoforms

`preprocess_annotation` takes per-sample TPMs (Salmon, Kallisto, FLAIR, or Mandalorian) plus the
GTF/GFF, computes max/min/median/mean TPM per transcript, and filters lowly-expressed isoforms (all
TPM > 0 kept by default). `--top_n N` keeps only the N most-expressed isoforms per gene. `--tss_window`
/ `--tes_window` emit BED windows for dCas targeting (feed them to `generate_guides -k` and
`rank_guides -t`).

```bash
crisprware preprocess_annotation -g chr19_ucsc_mm39.ncbiRefSeq.gtf \
  -t quant1.sf quant2.sf quant3.sf --type infer \
  --median 5 --top_n 10 --top_n_column median \
  --model consensus metagene shortest longest \
  --tss_window 300 300 --tes_window 300 300
```

```{important}
The GTF and the TPM files must share transcript IDs, and don't mix quantifiers (e.g. Salmon + Kallisto)
in one run.
```

## Ribo-seq: target translated ORFs

Convert ORF calls into a GTF with CDS entries, then run the pipeline as usual. For
[RiboTISH](https://github.com/zhpn1024/ribotish) run `ribotish predict` with `--inframecount --blocks
--aaseq` and pass the same GTF; [Price](https://github.com/erhard-lab/price) works with defaults.

```bash
gtf_from_ribotish.py -r ribotish_predict.tsv -i annotation.gtf -o orfs.gtf   # + filters
gtf_from_price.py    -i price.tsv             -g annotation.gtf -o orfs.gtf   # + filters
```

Both have filter flags (amino-acid length, in-frame count, TIS/frame/Fisher q-values, start codon,
...); see `gtf_from_ribotish.py -h` / `gtf_from_price.py -h`. For other ORF callers, open a GitHub issue.
