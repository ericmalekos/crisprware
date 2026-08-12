# Index and guides

Two preparatory builds: the off-target **index** every guide is scanned against, and the **guide set**
itself. Both differ by nuclease: AsCas12a uses the strong `TTTV` PAM for its guides and its off-target
index; enCas12a uses its broader high-activity `fa>0.6` on-target PAM set for guides and the even
broader `fa>0.3` set for the off-target index.

## 1. Off-target index

The index PAM depends on the enzyme you will build a track for.

::::{tab-set}

:::{tab-item} AsCas12a (default)
```bash
crisprware index_genome \
  -f $PROJ/hg38_primary_chr.fna \
  --pam TTTV \
  -l 23 \
  --pam_5_prime \
  --bin_width 14 \
  -o $PROJ
```
`--pam TTTV` is AsCas12a's off-target search space. AsCas12a does **not** cleave `TTTT`, so a `TTTV`
index is exact -- there are no weak-PAM off-targets to miss, and no genome PAM lookup is needed at
scoring time.
:::

:::{tab-item} enCas12a
```bash
crisprware index_genome \
  -f $PROJ/hg38_primary_chr.fna \
  --pam_list_file encas12a_offtarget_fa0.3.pam \
  -l 23 \
  --pam_5_prime \
  --bin_width 14 \
  -o $PROJ/encas12a
```
`--pam_list_file` unions an explicit PAM set (not a single IUPAC motif) into the index. Use the **fa>0.3**
set: the 30 four-mer PAMs enAsCas12a cleaves at >0.3 relative activity. It ships as
`encas12a_offtarget_fa0.3.pam` in the cas12a-screen-ot data and is, equivalently, the `in_fa0.3 == 1`
rows of the shipped `crisprware/data/encas12a_offtarget_weights.tsv`. enCas12a cleaves these weak PAMs,
so its off-target search must include them.
:::

::::

| Flag | Meaning |
|------|---------|
| `-l 23` | 23-nt protospacer length (Cas12a). |
| `--pam_5_prime` | PAM sits 5' of the spacer (Cas12a geometry). |
| `--bin_width 14` | index bin width; 14 is the measured throughput optimum that still fits a 24 GB GPU. |

Output: a directory `hg38_primary_chr_crisprots/` holding a **single** crispr-ots index file,
`hg38_primary_chr_crisprots.crot` (a few GB for hg38; the broad `fa>0.3` enCas12a index is larger than
the `TTTV` AsCas12a index). At scoring time you pass this as `-i` **without** the `.crot` extension, i.e.
the prefix `$PROJ/hg38_primary_chr_crisprots/hg38_primary_chr_crisprots`. Build is single-threaded-ish
and memory-light (~19 GB RAM, ~10-22 min for a human genome).

```{note}
`index_genome` names the output directory after the **input FASTA** (`<fasta_basename>_crisprots`),
independent of the PAM. Building the AsCas12a (`TTTV`) and enCas12a (`fa>0.3`) indexes from the same
FASTA into the same `-o` directory would therefore collide; build them into different `-o` directories
(above, the enCas12a index goes to `$PROJ/encas12a/...`).
```

## 2. Enumerate on-target protospacers

The guide PAM set is the nuclease's on-target search space, so it differs by enzyme.

::::{tab-set}

:::{tab-item} AsCas12a (default)
```bash
crisprware generate_guides \
  -f $PROJ/hg38_primary_chr.fna \
  --pam TTTV \
  --pam_5_prime -5 19 -3 23 -l 23 -w 4 3 \
  --threads 16 \
  -o $PROJ
```
`--pam TTTV` is the strong canonical AsCas12a PAM (`V = A/C/G`, excludes the weak `TTTT`).
:::

:::{tab-item} enCas12a
```bash
crisprware generate_guides \
  -f $PROJ/hg38_primary_chr.fna \
  --pam_list_file encas12a_ontarget.pam \
  --pam_5_prime -5 19 -3 23 -l 23 -w 4 3 \
  --threads 16 \
  -o $PROJ/encas12a
```
enCas12a targets a broader set of high-activity PAMs, so its guides use the **fa>0.6** on-target set:
the 15 four-mer PAMs enAsCas12a cleaves at >0.6 relative activity (`encas12a_ontarget.pam` in the
cas12a-screen-ot data; a superset of `TTTV`). Note the two enCas12a sets differ -- **fa>0.6** guides
but the broader **fa>0.3** off-target index above.
:::

::::

| Flag | Meaning |
|------|---------|
| `--pam_5_prime`, `-5 19 -3 23`, `-l 23` | Cas12a geometry: 5' PAM, 23-nt protospacer. |
| `-w 4 3` | context window: 4 nt up + 4 PAM + 23 spacer + 3 down = a **34-nt** window, exactly what the DeepCpf1 / EnPAM-GB / seq-DeepCpf1 models expect. |
| `--threads 16` | window-parallel engine; ~8x faster genome-wide. |

Output: `<fasta_basename>_gRNA/<fasta_basename>_gRNA.bed`, one row per guide. Build each enzyme's guides
into a separate `-o` directory, like the index, to avoid the FASTA-basename directory collision.
