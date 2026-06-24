# index_genome

Build the **crispr-ots** off-target index that `score_guides` scans against (fast, exact-mismatch).

```bash
crisprware index_genome -f tests/test_data/ce11/chrIII_sequence.fasta -p NGG -l 20
```

```text
chrIII_sequence_crisprots/
  chrIII_sequence_crisprots.crot      # the crispr-ots mmap off-target database
```

The index is PAM- and length-aware, so `-p/--pam` and `-l/--protospacer_length` are **required** and
must match the guides you will score (`NGG`/`20` for SpCas9, `TTTN`/`23 --pam_5_prime` for Cas12a).

## Options you reach for

- `--pam_5_prime` -- PAM is 5' of the protospacer (Cas12a).
- `--bin_width` -- crispr-ots index bin size (memory/throughput tradeoff).
- `--locations_to_keep <bed/gtf> [--feature CDS] [-w up down]` -- index only these regions (plus an
  optional flanking window) instead of the whole FASTA, for tighter, faster specificity.

```{tip}
For Cas12a use `-p TTTN -l 23 --pam_5_prime` (index `TTTN` so weak-PAM off-targets are still found).
Guidescan2 ships prebuilt indices for some genomes at <https://guidescan.com/downloads>.
```
