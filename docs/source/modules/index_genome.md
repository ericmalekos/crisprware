# index_genome

Build the **crispr-ots** off-target index that `score_guides` scans against (fast, exact-mismatch).

```bash
crisprware index_genome -f tests/test_data/ce11/chrIII_sequence.fasta -p NGG -l 20
```

```text
chrIII_sequence_crisprots/
  chrIII_sequence_crisprots.crot      # the crispr-ots mmap off-target database
```

The index is PAM- and length-aware, so `-l/--protospacer_length` plus **exactly one** of `-p/--pam`
or `--pam_list_file` are **required**, and must match the guides you will score (`NGG`/`20` for SpCas9,
`TTTN`/`23 --pam_5_prime` for Cas12a).

## Options you reach for

- `--pam_5_prime` -- PAM is 5' of the protospacer (Cas12a).
- `--pam_list_file <file>` -- index every PAM in a file (whitespace/comma/newline separated, `#`
  comments) instead of a single `-p/--pam` IUPAC pattern; use for PAM sets that are not one IUPAC
  motif (e.g. the broad enCas12a fa>0.3 off-target set). Mutually exclusive with `-p/--pam`.
- `--bin_width` -- crispr-ots index bin size (memory/throughput tradeoff).
- `--locations_to_keep <bed/gtf> [--feature CDS] [-w up down]` -- index only these regions (plus an
  optional flanking window) instead of the whole FASTA, for tighter, faster specificity.

```{tip}
For Cas12a use `-l 23 --pam_5_prime`. A plain `-p TTTN` indexes all TTT[ACGT]; the Cas12a **browser
tracks** need enzyme-specific indexes instead -- `-p TTTV` for AsCas12a and the broad fa>0.3 PAM set
(via `--pam_list_file`) for enCas12a (see the [Cas12a UCSC track pipeline](../ucsc_cas12a/index.md)).
Guidescan2 ships prebuilt indices for some genomes at <https://guidescan.com/downloads>.
```
