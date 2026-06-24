# Index and guides

Two preparatory builds: the off-target **index** every guide is scanned against, and the **guide set**
itself.

## 1. Off-target index (PAM `TTTN`)

```bash
crisprware index_genome \
  -f $PROJ/hg38_primary_chr.fna \
  --pam TTTN \
  -l 23 \
  --pam_5_prime \
  --bin_width 14 \
  -o $PROJ
```

| Flag | Meaning |
|------|---------|
| `--pam TTTN` | index PAM. `TTTN`, **not** `TTTV`, so weak-PAM (`TTTT`) off-targets are enumerable. |
| `-l 23` | 23-nt protospacer length (Cas12a). |
| `--pam_5_prime` | PAM sits 5' of the spacer (Cas12a geometry). |
| `--bin_width 14` | index bin width; 14 is the measured throughput optimum that still fits a 24 GB GPU. |

Output: `hg38_primary_chr_crisprots/` (the index, ~7 GB for hg38). Build is single-threaded-ish and
memory-light (~19 GB RAM, ~9-22 min for a human genome).

## 2. Enumerate `TTTV` protospacers

```bash
crisprware generate_guides \
  -f $PROJ/hg38_primary_chr.fna \
  --pam TTTV \
  --pam_5_prime \
  -5 19 -3 23 \
  -l 23 \
  -w 4 3 \
  --threads 16 \
  -o $PROJ
```

| Flag | Meaning |
|------|---------|
| `--pam TTTV` | guide PAM, the strong canonical set (`V = A/C/G`, excludes the weak `TTTT`). |
| `-5 19 -3 23`, `-l 23` | protospacer extent relative to the PAM. |
| `-w 4 3` | context window: 4 nt up + 4 PAM + 23 spacer + 3 down = a **34-nt** window, exactly what the DeepCpf1 / EnPAM-GB / enseq-DeepCpf1 models expect. |
| `--threads 16` | window-parallel engine; ~8x faster genome-wide. |

Output: `hg38_primary_chr_gRNA/hg38_primary_chr_gRNA.bed`, one row per guide.

```text
guides: 135,867,611
```

```{tip}
The index is `TTTN` but the guides are `TTTV`. That is deliberate (see the
[overview](00_overview.md#the-off-target-asymmetry-that-drives-the-index)): we score strong-PAM guides
against an index that can still find their weak-PAM off-targets.
```
