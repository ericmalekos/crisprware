# Command reference

Full `--help` for every `crisprware` subcommand. The [module pages](modules/index.md) cover the options you reach for most; this page is the complete listing.

## preprocess_annotation

```text
usage: crisprware preprocess_annotation [-h] -g GTF [-t [TPM_FILES ...]]
                                        [-f {salmon,kallisto,flair,mandalorian,infer}]
                                        [--mean MEAN] [--median MEDIAN]
                                        [--min MIN] [--max MAX] [-n TOP_N]
                                        [-c {median,mean,min,max}]
                                        [-m [{metagene,consensus,longest,shortest} ...]]
                                        [-w TSS_WINDOW TSS_WINDOW]
                                        [-e TES_WINDOW TES_WINDOW]
                                        [-x TX_TO_GENE] [--strip_tx_id]
                                        [-o OUTPUT_DIRECTORY]

Preprocess GTF/GFF annotations with optional RNA-seq filtering.

options:
  -h, --help            show this help message and exit
  -g GTF, --gtf GTF     GTF/GFF file to use for isoform filtering.
  -t [TPM_FILES ...], --tpm_files [TPM_FILES ...]
                        A list of one or more isoform quantification files
                        produced by Salmon, Kallisto or FLAIR (FLAIR outputs
                        counts, not TPMs). The first column should contain
                        only the transcript_id and should exactly match the
                        transcript_ids in --gtf. All transcript_ids in each
                        TPM file must be common across all files and must be
                        found in the GTF file.
  -f {salmon,kallisto,flair,mandalorian,infer}, --type {salmon,kallisto,flair,mandalorian,infer}
                        Specify TPM input type. 'infer' guesses the input type
                        based on the header line. [default: "infer"].
  --mean MEAN           For a given isoform, the mean tpm/count across samples
                        must be at least this to be considered, else discard
                        isoform. [default: 0.0]
  --median MEDIAN       For a given isoform, the median tpm/count across
                        samples must be at least this to be considered, else
                        discard isoform. [default: 0.0]
  --min MIN             For a given isoform, each sample must have at least
                        this tpm/count to be considered, else discard isoform.
                        [default: 0.0]
  --max MAX             For a given isoform, at least one sample must have at
                        least this tpm/count to be considered, else discard
                        isoform. [default: 0.0]
  -n TOP_N, --top_n TOP_N
                        For a given gene, rank all isoforms by median_tpm,
                        keep the top_n ranked isoforms and discard the rest.
                        '-1' to keep all isoforms. [default: -1]
  -c {median,mean,min,max}, --top_n_column {median,mean,min,max}
                        The metric by which to rank and filter top isoforms.
                        Used with '-n' to select expressed transcripts.
                        [default: median]
  -m [{metagene,consensus,longest,shortest} ...], --model [{metagene,consensus,longest,shortest} ...]
                        Whether to output 'metagene', 'consensus', 'longest',
                        'shortest' model. 'longest' and 'shortest' select, for
                        a given gene, the transcript with the longest or
                        shortest CDS, for now noncoding genes are ignored.
                        Output is always after tpm filtering has been applied.
                        Multiple entries are allowed e.g. --model metagene
                        consensus longest [default: None]
  -w TSS_WINDOW TSS_WINDOW, --tss_window TSS_WINDOW TSS_WINDOW
                        Pass two, space-separated, integers to specifiy the bp
                        window around the TSS as '<upstream>' '<downstream>'.
                        Strand-orientation is inferred, i.e. '<upstream>' will
                        be in the 5' direction of the TSS and <downstream> in
                        the 3' direction. e.g. --tss_window 250 150. [default:
                        None]
  -e TES_WINDOW TES_WINDOW, --tes_window TES_WINDOW TES_WINDOW
                        Pass two, space-separated, integers to specifiy the bp
                        window around the transcription end site, TES, as
                        '<upstream>' '<downstream>'. Strand-orientation is
                        inferred, i.e. '<upstream>' will be in the 5'
                        direction of the TES and <downstream> in the 3'
                        direction. e.g. --tss_window 0 150. [default: None]
  -x TX_TO_GENE, --tx_to_gene TX_TO_GENE
                        A TSV with transcript IDs in the first column and Gene
                        IDs in the second. The transcript IDs must match the
                        first column entries of the --quant_files. If this is
                        not provided it will be deduced from the GTF/GFF3 and
                        saved as './tmp/tx2gene.tsv'.
  --strip_tx_id         Set this flag if there are transcript IDs in the
                        quantification files but not in the GTF/GFF3.
                        [default: False]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Path to output. [default: current directory]
```

## index_genome

```text
usage: crisprware index_genome [-h] -f FASTA [-k [LOCATIONS_TO_KEEP ...]]
                               [--feature FEATURE]
                               [-w CONTEXT_WINDOW CONTEXT_WINDOW] -p PAM -l
                               PROTOSPACER_LENGTH [--pam_5_prime]
                               [--bin_width BIN_WIDTH] [-o OUTPUT_DIRECTORY]

Build a crispr-ots off-target index (PAM, protospacer length, and PAM
orientation define the enzyme).

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        FASTA file to use as a reference for index creation.
  -k [LOCATIONS_TO_KEEP ...], --locations_to_keep [LOCATIONS_TO_KEEP ...]
                        List of BED/GTF files with coordinates to use for
                        index creation. These locations will be used for off-
                        target scoring. If multiple files are passed,
                        coordinates will be merged with a union operation.
                        Leave empty to use entire fasta.
  --feature FEATURE     For any GTF/GFF in '--locations_to_keep', only this
                        feature will be used for determining appropriate
                        sgRNA. The feature should match an entry in the third
                        column of the GTF/GFF. [default: 'transcript']
  -w CONTEXT_WINDOW CONTEXT_WINDOW, --context_window CONTEXT_WINDOW CONTEXT_WINDOW
                        Pass two, space-separated, integers to specifiy the
                        nucleotide window around the --locations_to_keep
                        '<upstream>' '<downstream>'. This can be used to
                        expand the window around the final intervals e.g. '-w
                        1000 1500' expands chr1 2000 3500 -> chr1 1000 5000
                        Good for CRISPRi/a [default: 20 20]
  -p PAM, --pam PAM     PAM motif as an IUPAC string, e.g. NGG (SpCas9), NAG,
                        TTTV or TTTN (Cas12a). With --protospacer_length and
                        --pam_5_prime this fully defines the indexed off-
                        target sites.
  -l PROTOSPACER_LENGTH, --protospacer_length PROTOSPACER_LENGTH
                        Protospacer (spacer) length in bases, e.g. 20 for
                        SpCas9 or 23 for Cas12a.
  --pam_5_prime         Set if the PAM is 5' of the protospacer
                        (Cas12a-style). Omit for a 3' PAM (SpCas9-style).
                        [default: False]
  --bin_width BIN_WIDTH
                        crispr-ots bin-prefix width (1-15); higher = faster
                        GPU scan but a larger index. Sweet spots: ~9-10 SpCas9
                        full-human, 14 for Cas12a full-human on a 24 GB GPU.
                        [default: engine default]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Path to output. [default: current directory]
```

## generate_guides

```text
usage: crisprware generate_guides [-h] -f FASTA [-p PAM] [-l SGRNA_LENGTH]
                                  [-w CONTEXT_WINDOW CONTEXT_WINDOW]
                                  [-5 ACTIVE_SITE_OFFSET_5]
                                  [-3 ACTIVE_SITE_OFFSET_3]
                                  [-k [LOCATIONS_TO_KEEP ...]]
                                  [--feature FEATURE]
                                  [--join_operation {merge,intersect}]
                                  [--locations_to_discard [LOCATIONS_TO_DISCARD ...]]
                                  [--prefix PREFIX]
                                  [--gc_range GC_RANGE GC_RANGE]
                                  [--discard_poly_T] [--discard_poly_G]
                                  [--restriction_patterns [RESTRICTION_PATTERNS ...]]
                                  [--flank_5 FLANK_5] [--flank_3 FLANK_3]
                                  [--min_chr_length MIN_CHR_LENGTH]
                                  [--pam_5_prime] [--coords_as_active_site]
                                  [-t THREADS] [--chunk_size CHUNK_SIZE]
                                  [-o OUTPUT_DIRECTORY]

Generate sgRNA sequences matching specified PAM.

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        FASTA file to use as a reference for sgRNA generation.
  -p PAM, --pam PAM     Protospacer adjacent motif to match. All IUPAC
                        ambiguity codes are accepted as well as standard ATCG.
                        [default: NGG]
  -l SGRNA_LENGTH, --sgRNA_length SGRNA_LENGTH
                        Length of sgRNA to generate. [default: 20]
  -w CONTEXT_WINDOW CONTEXT_WINDOW, --context_window CONTEXT_WINDOW CONTEXT_WINDOW
                        Pass two, space-separated, integers to specifiy the
                        nucleotide window around the sgRNA as '<upstream>'
                        '<downstream>'. This can be used for downstream
                        scoring, For Ruleset3 use -w 4 6 to obtain an
                        appropriate score context. [default: 4 6]
  -5 ACTIVE_SITE_OFFSET_5, --active_site_offset_5 ACTIVE_SITE_OFFSET_5
                        Where cut occurs relative to PAM 5' end. To avoid
                        error, use '=' sign when passing a negative number,
                        e.g. --active_site_offset_5=-1 [default: -4]
  -3 ACTIVE_SITE_OFFSET_3, --active_site_offset_3 ACTIVE_SITE_OFFSET_3
                        Where cut occurs relative to PAM 5' end. [default: -2]
                        To avoid error, use '=' sign when passing a negative
                        number, e.g. --active_site_offset_3=-3 [default: -4]
  -k [LOCATIONS_TO_KEEP ...], --locations_to_keep [LOCATIONS_TO_KEEP ...]
                        List of BED/GTF files with coordinates in which the
                        sgRNA desired. If the sgRNA cutsite does not intersect
                        coordinates in these files they are discarded. Leave
                        blank to keep all sgRNA. e.g. atac_peak.bed genes.gtf
  --feature FEATURE     For any GTF/GFF in '--locations_to_keep', only this
                        feature will be used for determining appropriate
                        sgRNA. The feature should match an entry in the third
                        column of the GTF/GFF. [default: 'exon']
  --join_operation {merge,intersect}
                        How to treat '--locations_to_keep' if multiple files
                        are passed. Either 'merge' or 'intersect' can be used
                        and work as described in Bedtools. If 'merge', sgRNA
                        will be kept if its cutsite intersects an entry in ANY
                        of the files, if 'intersect' the cutsite must
                        intersect an entry in EACH file. [default:
                        'intersect']
  --locations_to_discard [LOCATIONS_TO_DISCARD ...]
                        List of BED/GTF files with coordinates where sgRNA
                        should not target. If the sgRNA cutsite intersects
                        coordinates in these files the sgRNA is discarded.
                        Leave blank to keep all sgRNA. e.g. TSS.bed
                        coding_genes.gtf
  --prefix PREFIX       Prefix to use for sgRNA identifiers. [default: None]
  --gc_range GC_RANGE GC_RANGE
                        Pass two, space-separated, integers to specifiy the
                        percentile range of GC content e.g. '--gc_range 25
                        75'. [default: 0 100]
  --discard_poly_T      Whether to discard polyT (>TTT) sgRNA. Recommend True
                        for PolIII promoters [default: False]
  --discard_poly_G      Whether to discard polyT (>GGGG) sgRNA. [default:
                        False]
  --restriction_patterns [RESTRICTION_PATTERNS ...]
                        Reject sgRNA with these restriction patterns. Also
                        checks 5'flank+sgRNA+3'flank, and reverse complement,
                        if provided. For multiple values, separate by space.
                        e.g. GCGGCCGC TCTAGA CACCTGC
  --flank_5 FLANK_5     include the 5' context of the lentivirus vector. Used
                        in conjunction with --restriction_patterns to remove
                        incompatible sgRNA
  --flank_3 FLANK_3     include the 3' context of the lentivirus vector. Used
                        in conjunction with --restriction_patterns to remove
                        incompatible sgRNA
  --min_chr_length MIN_CHR_LENGTH
                        Minimum chromosome length to consider for sgRNA
                        generation. [default: 20]
  --pam_5_prime         If the PAM is positioned 5' to the protospacer set
                        this flag, e.g. for Cas12a sgRNAs [default: False]
  --coords_as_active_site
                        Whether to output bed coordinates at the active site
                        rather than the coordinates of the entire protospacer.
                        For purposes of keeping or discarding sgRNAs, overlap
                        with the active site coordinates will be used
                        regardless [default: True]
  -t THREADS, --threads THREADS
                        Number of worker processes. Parallelism scales with
                        the number of genome windows (see --chunk_size), not
                        the PAM count, so this can usefully exceed it.
                        [default: 4]
  --chunk_size CHUNK_SIZE
                        Genome window size in bp for parallel scanning;
                        smaller = more/finer tasks (better load balance,
                        slightly more overhead). [default: 5,000,000]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Path to output. [default: current directory]
```

## score_guides

```text
usage: crisprware score_guides [-h] -b GRNA_BED [-i [GUIDESCAN2_INDICES ...]]
                               [--tracr {Hsu2013,Chen2013,both}] [-t THREADS]
                               [--threshold THRESHOLD]
                               [--mismatches MISMATCHES]
                               [--rna_bulges RNA_BULGES]
                               [--dna_bulges DNA_BULGES]
                               [--mode {succinct,complete}]
                               [--alt_pams [ALT_PAMS ...]] [-d] [--skip_rs3]
                               [--skip_gs2] [--min_rs3 MIN_RS3]
                               [--cas12a_scorer {none,enpam_gb,deepcpf1,enseq_deepcpf1,seq_deepcpf1variants,both} [{none,enpam_gb,deepcpf1,enseq_deepcpf1,seq_deepcpf1variants,both} ...]]
                               [--cas12a_variant CAS12A_VARIANT]
                               [--min_deepcpf1 MIN_DEEPCPF1]
                               [--min_enpam_gb MIN_ENPAM_GB]
                               [--min_enseq_deepcpf1 MIN_ENSEQ_DEEPCPF1]
                               [--cas9_scorer {none,deepspcas9,deephf_wt_u6,deephf_esp,deephf_hf} [{none,deepspcas9,deephf_wt_u6,deephf_esp,deephf_hf} ...]]
                               [--min_deepspcas9 MIN_DEEPSPCAS9]
                               [--min_deephf MIN_DEEPHF]
                               [--chunk_size CHUNK_SIZE] [-k]
                               [-o OUTPUT_DIRECTORY] [--ucscgb UCSCGB]
                               [--chrom_sizes CHROM_SIZES]
                               [--crispr_ots_bin CRISPR_OTS_BIN]
                               [--ucscgb_scanner {gpu,cpu}]
                               [--ucscgb_cfd_threshold UCSCGB_CFD_THRESHOLD]
                               [--ucscgb_list_cap UCSCGB_LIST_CAP]
                               [--ucscgb_blank_threshold UCSCGB_BLANK_THRESHOLD]
                               [--ucscgb_reuse_offtargets]
                               [--ucscgb_2xnls | --no-ucscgb_2xnls]
                               [--ucscgb_chunk_max UCSCGB_CHUNK_MAX]
                               [--ucscgb_gpus UCSCGB_GPUS]
                               [--ucscgb_keep_chunks]

Score guides with RS3 cleavage and crispr-ots/Guidescan2 off-target
specificity.

options:
  -h, --help            show this help message and exit
  -b GRNA_BED, --grna_bed GRNA_BED
                        grnas.bed ouput of generate_guides.
  -i [GUIDESCAN2_INDICES ...], --guidescan2_indices [GUIDESCAN2_INDICES ...]
                        One or more, space-separate Guidescan2 indices. A
                        specificity score will be calculated against each
                        index separately.
  --tracr {Hsu2013,Chen2013,both}
                        TracrRNA version for cleavage scoring. Either
                        'Hsu2013' or 'Chen2013' or 'both', see
                        https://github.com/gpp-rnd/rs3 for details.
  -t THREADS, --threads THREADS
                        Number of threads [default: 8]
  --threshold THRESHOLD
                        Threshold for Guidescan2 off-target hits. If off-
                        targets are found this distance away the gRNA will be
                        discarded, i.e. set to 2 to discard any guides with a
                        0, 1 or 2 mismatches from another PAM adjacent
                        sequence. --threshold=-1 to retain all guides
                        [default: 2]
  --mismatches MISMATCHES
                        Number of mismatches for Guidescan2 off-target scoring
                        [default: 3]
  --rna_bulges RNA_BULGES
                        RNA bulges for Guidescan2 off-target scoring [default:
                        0]
  --dna_bulges DNA_BULGES
                        DNA bulges for Guidescan2 off-target scoring [default:
                        0]
  --mode {succinct,complete}
                        Whether Guidescan2 temporary output should be succinct
                        or complete mode [default: 0]
  --alt_pams [ALT_PAMS ...]
                        One or more, space-separate alternative pams for off-
                        target consideration. e.g. NAG
  -d, --drop_duplicates
                        Drop exact duplicate gRNAs before scoring to save
                        time. Set flag to retain duplicates. [default: True]
  --skip_rs3            Set flag to skip RS3 scoring [default: False]
  --skip_gs2            Set flag to skip Guidescan2 scoring [default: False]
  --min_rs3 MIN_RS3     Minimum cleavage RS3 score. RS3 cleavage scores are
                        formatted as z-scores, so this is interpreted as a
                        standard deviation cutoff. Functionality also
                        available in rank_guides.py. Applying at this stage
                        can increase speed by filtering before off-target
                        scoring. [default: None]
  --cas12a_scorer {none,enpam_gb,deepcpf1,enseq_deepcpf1,seq_deepcpf1variants,both} [{none,enpam_gb,deepcpf1,enseq_deepcpf1,seq_deepcpf1variants,both} ...]
                        One or more Cas12a on-target scorers (space-
                        separated). enpam_gb for en(As)Cas12a; deepcpf1 for
                        wildtype AsCas12a/LbCas12a (Kim 2018); enseq_deepcpf1
                        for wildtype AsCas12a (Chen 2025, modern);
                        seq_deepcpf1variants for variant-specific scoring
                        (requires --cas12a_variant). 'both' is the legacy
                        alias for 'enpam_gb deepcpf1'. Multiple values may be
                        combined (e.g. --cas12a_scorer enpam_gb deepcpf1
                        enseq_deepcpf1). Mutually exclusive with --tracr (RS3
                        is SpCas9-only); implies --skip_rs3. [default: none]
  --cas12a_variant CAS12A_VARIANT
                        Cas12a variant name (e.g. AsCas12a_Ultra,
                        enAsCas12a-HF1, LbCas12a, HyperFi-AsCas12a). Required
                        when --cas12a_scorer is seq_deepcpf1variants. See
                        crisprware.scorers.seq_deepcpf1variants for the full
                        23-variant list.
  --min_deepcpf1 MIN_DEEPCPF1
                        Minimum DeepCpf1 score (raw regression, ~[0, 100]).
                        Applied after scoring; analogous to --min_rs3.
                        [default: None]
  --min_enpam_gb MIN_ENPAM_GB
                        Minimum enPAM+GB score (probability-like, [0, 1]).
                        Applied after scoring; analogous to --min_rs3.
                        [default: None]
  --min_enseq_deepcpf1 MIN_ENSEQ_DEEPCPF1
                        Minimum enseq-DeepCpf1 / seq-DeepCpf1variants score
                        (probability, [0, 1]). Applied after scoring.
                        [default: None]
  --cas9_scorer {none,deepspcas9,deephf_wt_u6,deephf_esp,deephf_hf} [{none,deepspcas9,deephf_wt_u6,deephf_esp,deephf_hf} ...]
                        One or more additional SpCas9 on-target scorers
                        (space-separated) to run alongside RS3. deepspcas9:
                        Kim 2019 inception-CNN (Sci Adv); 30-nt context (4 +
                        20 protospacer + 3 PAM + 3 downstream); unbounded
                        regression. deephf_*: Wang 2019 BiLSTM (Nat Commun)
                        for three Cas9 variants -- wildtype SpCas9 (wt_u6),
                        eSpCas9 (esp), SpCas9-HF1 (hf); 23-nt protospacer+PAM
                        input; output in [0, 1]. Multiple values may be
                        combined (e.g. --cas9_scorer deepspcas9 deephf_wt_u6
                        deephf_esp deephf_hf). Runs in parallel with --tracr
                        (no mutex). [default: none]
  --min_deepspcas9 MIN_DEEPSPCAS9
                        Minimum DeepSpCas9 score (unbounded regression, ~[0,
                        100]). Applied after scoring; analogous to --min_rs3.
                        [default: None]
  --min_deephf MIN_DEEPHF
                        Minimum DeepHF score (probability, [0, 1]). Applied to
                        whichever deephf_* variant is selected. [default:
                        None]
  --chunk_size CHUNK_SIZE
                        Number of gRNAs to hold in memory for cleavage scoring
                        and off-target filtering. Reduce if memory
                        constrained. Increasing may improve runtime [default:
                        100000]
  -k, --keep_tmp        Set flag to keep temporary Guidescan2 output [default:
                        False]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Path to output. [default: current directory]

UCSC browser track (Cas12a):
  --ucscgb UCSCGB       Alternative output path: write a UCSC Genome Browser
                        Cas12a track (cas12a.bb, crisprDetails.tab.gz + .gzi,
                        cas12aTargets.as) to this directory. Runs the selected
                        --cas12a_scorer on-target models plus a streaming
                        crispr-ots off-target pass (--output-mode both,
                        skipping off-targets for perfect-match guides) against
                        the first -i index. Implies Cas12a mode and retains
                        duplicate guides. [default: None]
  --chrom_sizes CHROM_SIZES
                        chrom.sizes file (chrom<TAB>size) for bedToBigBed when
                        --ucscgb is set.
  --crispr_ots_bin CRISPR_OTS_BIN
                        Path to the crispr-ots binary used for the --ucscgb
                        streaming pass (must support --output-mode/--scanner;
                        default: 'crispr-ots' on PATH).
  --ucscgb_scanner {gpu,cpu}
                        Scanner for the --ucscgb off-target pass [default:
                        gpu].
  --ucscgb_cfd_threshold UCSCGB_CFD_THRESHOLD
                        CFD floor for the off-target list written to
                        crisprDetails (counts are unaffected — they come from
                        the unfloored Mode-1 totals). [default: 0.023]
  --ucscgb_list_cap UCSCGB_LIST_CAP
                        Max off-targets listed per guide in crisprDetails (top
                        by score) [default: 100].
  --ucscgb_blank_threshold UCSCGB_BLANK_THRESHOLD
                        Guides with more than this many off-targets get an
                        empty list (counts kept) — the viewer shows 'Too many
                        off-targets'. [default: 2000]
  --ucscgb_reuse_offtargets
                        Reuse an existing off-target enumerate output
                        (offtargets.csv + .ot.tsv) in the --ucscgb dir instead
                        of re-running crispr-ots. Use when re-assembling a
                        track after an adapter-only change (the engine output
                        is unchanged), skipping the GPU pass. [default: off]
  --ucscgb_2xnls, --no-ucscgb_2xnls
                        Also score off-targets with the 2xNLS-Cas12a matrix
                        and emit its TTTV/TTTN specificity columns (a second,
                        aggregated off-target pass — the off-target list still
                        comes from the enCas12a pass). [default: on; use --no-
                        ucscgb_2xnls to skip]
  --ucscgb_chunk_max UCSCGB_CHUNK_MAX
                        If >0, build the track in chunks of at most this many
                        guides and merge the per-chunk tracks into one
                        (recomputing crisprDetails.tab byte offsets). Bounds
                        peak memory and disk, and enables multi-GPU fan-out
                        via --ucscgb_gpus. 0 = single pass over all guides.
                        [default: 0]
  --ucscgb_gpus UCSCGB_GPUS
                        Comma-separated CUDA device ids to fan
                        --ucscgb_chunk_max chunks across (e.g.
                        '0,1,2,3,4,5,6,7'); each concurrent chunk's off-target
                        pass is pinned to one device, so chunks never share a
                        GPU. Only used when --ucscgb_chunk_max>0. [default:
                        single device 0].
  --ucscgb_keep_chunks  Keep the per-chunk working dir (<ucscgb>/chunks/)
                        after merging, for inspection or resume. [default:
                        delete after a successful merge].
```

## rank_guides

```text
usage: crisprware rank_guides [-h] -k SCORED_GUIDES -t TARGETS
                              [--target_mode {gene,tx}] [-f FEATURE]
                              [-p PERCENTILE_RANGE PERCENTILE_RANGE]
                              [-n NUMBER_OF_GUIDES]
                              [--min_spacing MIN_SPACING] [--output_all]
                              [--plot_histogram] [-c [FILTERING_COLUMNS ...]]
                              [-m [MINIMUM_VALUES ...]]
                              [-r [RANKING_COLUMNS ...]]
                              [-w [COLUMN_WEIGHTS ...]] [--normalize_columns]
                              [-o OUTPUT_DIRECTORY]

Rank and select best guides based on scoring criteria.

options:
  -h, --help            show this help message and exit
  -k SCORED_GUIDES, --scored_guides SCORED_GUIDES
                        <score_guides_output>.bed output from score_guides.
  -t TARGETS, --targets TARGETS
                        BED/GTF/GFF used to select final guides per target.
                        For GTF/GFF, set --target_mode to either 'gene' or
                        'tx'. For BED, targets are each entry. Use '--
                        number_of_targets' to set the number of guides chosen
                        for each target.
  --target_mode {gene,tx}
                        If a GTF/GFF is used to select targets, gRNAs can be
                        grouped at either the 'tx' or 'gene' level e.g. '--
                        target_mode gene -n 10' chooses 10 guides per gene, '
                        --target_mode tx -n 10' chooses 10 per transcript
                        [default: gene].
  -f FEATURE, --feature FEATURE
                        If GTF/GFF passed, use this feature for processing
                        e.g. 'exon', 'CDS', '5UTR', etc. The feature appears
                        in the third column of the GTF/GFF [default: CDS].
  -p PERCENTILE_RANGE PERCENTILE_RANGE, --percentile_range PERCENTILE_RANGE PERCENTILE_RANGE
                        Allowable range of guide for each transcript and
                        feature set, e.g. '-p 60 80 -f exon' returns gRNAs in
                        the 60th to 80th percentile of exons for a given
                        transcript. Default setting returns guides anywhere in
                        the CDS for each transcript [default: 0 100]
  -n NUMBER_OF_GUIDES, --number_of_guides NUMBER_OF_GUIDES
                        Number of guides returned per target.'-1' to keep all
                        guides [default: -1]
  --min_spacing MIN_SPACING
                        The minimum nucleotide space between guides for a
                        given target. e.g. --min_spacing 10, requires guides
                        10 nts appart. 0 to allow overlapping guides.[default:
                        0]
  --output_all          Set flag to save gRNA-target TSVs at each stage of
                        filtering rather than just the end.[default: False]
  --plot_histogram      Set flag to plot a histogram of the distribution of
                        gRNAs per target after each filtering step. Sets '--
                        output_all' to True.[default: False]
  -c [FILTERING_COLUMNS ...], --filtering_columns [FILTERING_COLUMNS ...]
                        One or more space-separated column names used for
                        filtering. Uses raw values. e.g. '-c rs3_z_score
                        specificity_Hg38_index'.
  -m [MINIMUM_VALUES ...], --minimum_values [MINIMUM_VALUES ...]
                        A space-separated list of minimum values for each
                        column in passed by --ranking_columns. e.g. '-c
                        rs3_z_score specificity_Hg38_index -m "-1" 0.2'
                        Default is no minimum [default: None]
  -r [RANKING_COLUMNS ...], --ranking_columns [RANKING_COLUMNS ...]
                        One or more space-separated column names used for
                        guide ranking. e.g. '-r rs3_score_Hsu2013
                        rs3_score_Chen2013'.
  -w [COLUMN_WEIGHTS ...], --column_weights [COLUMN_WEIGHTS ...]
                        A space-separated list of weight values for each
                        column in passed by --ranking_columns. e.g. '-c
                        rs3_score specificity_Hg38_index -w 1 0' Default is
                        equal weighting for all ranking columns.
  --normalize_columns   Scale ranking column values to 0 to 1 [default: True]
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Path to output. [default: current directory]
```
