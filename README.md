# CRISPRware

CRISPRware is a comprehensive toolkit designed to preprocess NGS data and identify, score, and rank guide RNAs (gRNAs) for CRISPR experiments. It supports RNASeq, RiboSeq, ATACSeq, DNASESeq, ChIPSeq, and other genomic preprocessing techniques.

## Table of Contents
1. [Installation](#installation)
2. [Tutorials](#tutorials)
3. [Quickstart](#usage)
4. [Requirements](#requirements)
5. [Leveraging NGS data](#methods)
    - [RNASeq Guided Preprocessing](#rnaseq-guided-preprocessing)
    - [RiboSeq Guided Preprocessing](#riboseq-guided-preprocessing)
6. [Alternate PAMs and scoring methods](#altpams-and-scoring)
7. [References](#references)
8. [Full Commands](#commands)

## Installation

If you have not already, install one of the package managers [miniconda](https://docs.anaconda.com/miniconda/) or [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)

### Linux installation
With conda installed perform the following commands.
(If you installed micromamba, conda -> micromamba)

```
git clone https://github.com/ericmalekos/crisprware crisprware && cd crisprware

conda env create -f environment.yml && conda activate crisprware

pip install .
```

### MacOs installation and troubleshooting
Try running `git -h`, if you hit an error `xcrun: error: invalid developer path ...`, you may need to install the Command Line Tools package with `xcode-select --install`
With this complete, follow the same instructions as for Linux.

You may encounter an error with the `score_guides` module. In short you need to install a specific version of `libomp` for RS3 scoring, this can be accomplished with the command `brew install libomp@11.1.0`. You need to have [homebrew](https://brew.sh/) installed. 

### Docker 
Avoid local installation by pulling the latest docker image
and run commands:
```
docker pull ericmalekos/crisprware:latest

docker run crisprware preprocess_gtf -h
```

## Tutorials
These interactive notebooks demonstrate use of CRISPRware modules with text explanations and codeblocks. These are the best place to start for gaining understanding of the workflow and capabilities of this software. In each case the first block sets up the environment by pulling the latest version from Github.

### Full Tutorial: <a href="https://colab.research.google.com/github/ericmalekos/crisprware/blob/main/colab_notebooks/CRISPRware_Tutorial.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab" style="vertical-align: middle; margin-left: 5px;"/></a>
Covers all major CRISPRware functions.

### CRISPRware Rice Genome Tutorial: <a href="https://colab.research.google.com/github/ericmalekos/crisprware/blob/main/colab_notebooks/RiceGenome_Tutorial.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab" style="vertical-align: middle; margin-left: 5px;"/></a>
End to end processing of rice osa1_r7 genome and gene annotation.

### CRISPRware NGS applications: <a href="https://colab.research.google.com/github/ericmalekos/crisprware/blob/main/colab_notebooks/NGS_Tutorial.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab" style="vertical-align: middle; margin-left: 5px;"/></a>
Examples of retrieving public NGS data and applying to custom gRNA library design.

## Quickstart
### Input Requirements
- **FASTA File**
### Optional inputs
- **BED file**: A BED file can be provided to specify regions of interest within the genome. This file can help to limit the search space for gRNA identification.
- **GTF/GFF file**: A GTF or GFF file can be used to provide gene annotations. This information can be used to filter gRNAs based on specific genomic features such as exons or coding sequences .

### CRISPRware workflow
We demonstrate usage with ce11 chromosome III fasta and NCBI GTF, included in the tests/test_data/ce11 directory: 

Note the example off-target index is limited to chrIII, not the full ce11 genome

```
index_genome -f tests/test_data/ce11/chrIII_sequence.fasta
```

We can build gene models from NCBI GTF,

```
preprocess_annotation -g tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf \
-m metagene consensus longest shortest
```

Default settings generate NGG protospacer guides

```
generate_guides -f tests/test_data/ce11/chrIII_sequence.fasta \
-k tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf \
--feature CDS
```

Scoring will take ~5 minutes and uses 8 threads by default.
Change this with `--threads` <int>. `--tracr` is either `Chen2013`,`Hsu2013`, os `both`, see [RuleSet3](https://github.com/gpp-rnd/rs3) scoring for details

```
score_guides -b chrIII_sequence_gRNA/chrIII_sequence_gRNA.bed \
-i chrIII_sequence_gscan2/chrIII_sequence_gscan2 --tracr Chen2013 --threads 8
```

Ranking is done based on scoring columns  
`-c` is matched with `-m` order so this filters out `RS3_score_Chen2013 < 0`,  `specificity_gscan_index < 0.2`  
`-p 5 65`, `-f CDS` filters out gRNAs outside of the 5th-65th percentile of the CDS  
`--output_all` outputs TSV and histograms for each stage of filtering in addition to the final output.

```
rank_guides \
-k chrIII_sequence_scoredgRNA/chrIII_sequence_scoredgRNA.bed \
-t tests/test_data/ce11/chrIII_ce11.ncbiRefSeq.gtf \
-f CDS \
-c RS3_score_Chen2013 specificity_chrIII_sequence_gscan2 \
-m 0 0.2 \
-p 5 65 \
-r RS3_score_Chen2013 \
--output_all
```

## Requirements

Memory requirements may be substantial in both the `index_genome` and `score_guides` steps. Guidescan2 authors provide compiled indices for some model species in the [download section of their website](https://guidescan.com/downloads) which can be downloaded directly to avoid use of `index_genome`.  
For `score_guides` we provide a parameter `--chunk_size <n>` which can be used to decrease memory usage by processing `<n>` guides at a time instead of all at once. Default setting is 100000. Increasing this number will speed up processing time and memory requirements, decreasing will slow down processing time and decrease memory requirements.

## Leveraging NGS data

CRISPRware offers a series of modules to preprocess NGS data and determine suitable gRNAs for CRISPR applications.

### RNASeq Guided Preprocessing

The module `preprocess_annotation` takes processed RNASeq TPMs from Kallisto, Salmon, FLAIR, or Mandalorian from one or more samples along with the GTF/GFF gene annotation. All processed samples should be from the same quantification tool, don't mix Salmon and Kallisto files. If multiple samples are passed, max, min, median, and mean TPM values for each transcript are determined, and the user can supply minimum cut-offs for any combination of these to filter out lowly expressed isoforms. All detected isoforms (TPM > 0) are kept by default. The user can also set an integer flag `--top-n <n>` which will filter out all but the <n> most highly expressed isoform for each gene. So, `--top-n 1` will retain only the gene model of the most highly expressed isoform - according to median_tpm if multiple RNA seq files are passed. 
There are also `--tss_window` and `--tes_window` options, which produce BED for dCas target choices. User can use these GTFs/BEDs in the `generate_guides` step and the `rank_guides` step. 

```
preprocess_annotation -g test_data/chr19_ucsc_mm39.ncbiRefSeq.gtf \
-t quant1.sf quant2.sf quant3.sf \
--type infer \
--median 5 \
--top_n 10 \
--top_n_column median \
--model consensus metagene shortest longest \
--tss_window 300 300
--tes_window 300 300
```

*IMPORTANT: ensure the GTF and the TPM files have the same transcript IDs*


### RiboSeq Guided Preprocessing

A number of tools exists for calling translated ORFs from RiboSeq. In order to find gRNAs against these putative coding regions we can convert output from these programs into a GTF with annotated coding sequence (CDS) entries and run the pipeline normally. 

For ORFs called with [RiboTISH](https://github.com/zhpn1024/ribotish/tree/master) set these options in the `ribotish predict` command: `--inframecount`, `--blocks`, `--aaseq` and provide the same GTF that was passed to `ribotish`. Default settings should work for ORFs called with [Price](https://github.com/erhard-lab/price), but it does have fewer filtering options.

For other RiboSeq ORF callers raise a github issue and I will address it.

Full filtering options:

```
gtf_from_ribotish.py -h

options:
  -h, --help            show this help message and exit
  -r RIBOTISH, --ribotish RIBOTISH
                        Path to the Ribotish predict TSV file
  -i INPUT_GTF, --input_gtf INPUT_GTF
                        Path to the corresponding GTF file
  -o OUTPUT_GTF, --output_gtf OUTPUT_GTF
                        Path to output the new GTF file
  --min_aalen MIN_AALEN
                        Minimum amino acid length
  --min_inframecount MIN_INFRAMECOUNT
                        Minimum in-frame count
  --max_tisqvalue MAX_TISQVALUE
                        Maximum TIS Q-value
  --max_frameqvalue MAX_FRAMEQVALUE
                        Maximum Frame Q-value
  --max_fisherqvalue MAX_FISHERQVALUE
                        Maximum Fisher Q-value
  --select_based_on {AALen,InFrameCount,TISQvalue,FrameQvalue,FisherQvalue}
                        Column to select the best row for each Tid, TisType pair
  --genetype GENETYPE   GeneType to filter, must match a column entry
  --tistype TISTYPE     TisType to filter, must match a column entry


gtf_from_price.py -h

options:
  -h, --help            show this help message and exit
  -i INPUT_TSV, --input_tsv INPUT_TSV
                        Path to the input price TSV file
  -g INPUT_GTF, --input_gtf INPUT_GTF
                        Path to the input GTF file to be used as a reference
  -o OUTPUT_GTF, --output_gtf OUTPUT_GTF
                        Path to output the new GTF file
  -p MIN_P_VALUE, --min_p_value MIN_P_VALUE
                        Minimum p value for filtering
  --min_aalen MIN_AALEN
                        Minimum amino acid length
  --tis_type TIS_TYPE   Tis Type to filter
  --start_codon START_CODON
                        start codon to filter

```

## Alternate PAMs and scoring methods

Default `generate_guides` settings are equivalent to

```
generate_guides \
-f <fasta> \
--pam [-p] NGG
--sgRNA_length [-l] 20
--context_window [-w] 4 6
--active_site_offset_5 [-5] "-4"
--active_site_offset_3 [-3] "-4"

```

![plot](./images/NGG_cleavage.png)

All [IUPAC ambiguity codes](https://genome.ucsc.edu/goldenPath/help/iupac.html) are allowed and will be automically expanded, e.g. NGG -> AGG, TGG, CGG, GGG. Note that `context_window[0]` extends the sequence in the 5' direction, `context_window[1]` in the 3' direction. `active_site_offset`s are calculated relative to PAM-protospacer position, and should be passed in quotes if they are negative.


For Cas12A guide selection change `generate_guides` settings to

```
generate_guides \ 
-f <fasta> \
--pam TTTV --pam_5_prime -5 19 -3 23 -l 23 -w 8 2
```


![plot](./images/TTTV_cleavage.png)

Here the pam is 5-prime to the protospacer so `--pam_5_prime` flag is set and the length is increased 23. The window is resized for compatibility with DeepCpf1 and EnPAMGB scoring and final sequence should be 34 nts long.

### Additional scoring methods

For additional on-target scoring, including of Cas12A/Cpf1 guides, first install [crisprScore](https://github.com/crisprVerse/crisprScore) (recommendation: install in a new conda environment). Once installed the crisprScore.R script can be used to score guides.

```
conda activate <crisprscore env>

./scripts/crisprscore.R

	Usage: crisprscore.R <path_to_sgrna_bed_file> <method_number> <outputfile> [<additional settings> ... ] 
	Example: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 2 Azimuth_scored_sgRNAs.bed
	Methods:
	1:  RuleSet1 - SpCas9 (Length: 30)
	2:  Azimuth - SpCas9 (Length: 30)
	3:  DeepHF - SpCas9 (Length: 23)
	4:  Lindel - SpCas9 (Length: 65)
	5:  DeepCpf1 - AsCas12a (Length: 34)
	6:  EnPAMGB - enAsCas12a (Length: 34)
	7:  CRISPRscan - SpCas9 (Length: 35)
	8:  CasRx-RF - CasRx (Length: NA)
	9:  CRISPRai - SpCas9 (Length: 22)
	10: CRISPRater - SpCas9 (Length: 20)
	11: DeepSpCas9 - SpCas9 (Length: 30)
	12: RuleSet3 - SpCas9 (Length: 30)


	Additional optional argument for chunk size:
	--chunk-size: Specify the size of chunks for processing the dataframe (optional)
	Example with chunk size: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 2 Azimuth_scored_sgRNAs.bed --chunk-size 10000

	Additional settings for method 3 (DeepHF):
	enzyme: Specify the enzyme (options: 'WT', 'ESP', 'HF')
	promoter: Specify the promoter (options: 'U6', 'T7')
	Example: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 3 DeepHF_scored_sgRNAs.bed WT U6

	Additional setting for method 5 (DeepCpf1):
	--no-convertPAM: Specify whether non-canonical PAMs are converted to TTTC [default: TRUE]
	Example: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 5 DeepCpf1_scored_sgRNAs.bed --no-convertPAM

	Additional setting for method 12 (RuleSet3):
	--tracrRNA: Specify tracrRNA (options: 'Hsu2013', 'Chen2013')

	Example: crisprscore.R tests/test_data/chr19_GRCm39_sgRNA.bed 12 rs3_scored_sgRNAs.bed Chen2013
```

The output of this script should be passed to `score_guides` in order to properly format for `rank_guides`. Additional score columns will be added unless the user specifies `--skip_rs3` and/or `--skip_gs2`.

```
crisprscore.R 34_sgrnas.bed 6 EnPAMGB_sgRNAs.bed

score_guides -b EnPAMGB_sgRNAs.bed --skip_rs3 --skip_gs2
```


Guidescan2 is not compatible with PAMs 5' to protospacers, for off-target scoring in these cases I suggest [Flash Fry](https://github.com/mckennalab/FlashFry?tab=readme-ov-file)

## References
When using the CRISPRware track hub in your research, please cite:

<ul>
    <li><a href="https://www.biorxiv.org/content/10.1101/2024.06.18.599405v1">CRISPRware</a></li>
    <li><a href="https://www.biorxiv.org/content/10.1101/2022.05.02.490368v1">Guidescan2</a></li>
    <li><a href="https://www.nature.com/articles/s41467-022-34320-7">crispVerse</a></li>
</ul>

<p>And the cleavage-score method(s) you used:</p>

<ul>
    <li><a href="https://www.nature.com/articles/nbt.3026">Ruleset1</a></li>
    <li><a href="https://www.nature.com/articles/nbt.3437">Azimuth</a></li>
    <li><a href="https://www.nature.com/articles/s41467-022-33024-2">Ruleset3</a></li>
    <li><a href="https://academic.oup.com/nar/article/46/3/1375/4754467">CRISPRater</a></li>
    <li><a href="https://www.nature.com/articles/nmeth.3543">CRISPRscan</a></li>
    <li><a href="https://www.nature.com/articles/s41467-019-12281-8">DeepHF</a></li>
    <li><a href="https://www.science.org/doi/10.1126/sciadv.aax9249">DeepSpCas9</a></li>
    <li><a href="https://www.nature.com/articles/s41587-020-0600-6">EnPamGB</a></li>
    <li><a href="https://www.nature.com/articles/nbt.4061">DeepCpf1</a></li>
</ul>