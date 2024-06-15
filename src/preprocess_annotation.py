#!/usr/bin/env python3

'''

'''
import argparse
from utils.utility_functions import create_output_directory
from utils.gtf_bed_processing_functions import create_metagene_model,\
    create_constitutive_model,filter_gtf_by_transcript_ids, parse_gtf_for_cds_extremes,\
    extract_transcript_gene_relationship, gtf_to_tss_tes_bed
from utils.quantified_rna_functions import add_gene_ids_and_subset,\
    filter_dataframe,process_files

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Preprocess GTF for improved sgRNA selection."
    )

    parser.add_argument(
        "-g", "--gtf",
        type=str,
        help="GTF/GFF file to use for isoform filtering.",
        required=True
    )

    parser.add_argument(
        "-t", "--tpm_files",
        help="A list of one or more isoform quantification files produced by Salmon, \
        Kallisto or FLAIR (FLAIR outputs counts, not TPMs). The first column \
        should contain only the transcript_id and should exactly match the \
        transcript_ids in --gtf. All transcript_ids in each TPM file must be \
        common across all files and must be found in the GTF file.",
        type=str,
        default="",
        nargs='*'
    )

    parser.add_argument(
        "-f", "--type",
        type=str,
        choices=["salmon", "kallisto", "flair", "mandalorian", "infer"],
        help="Specify TPM input type. 'infer' guesses the input type based on the header line. [default: \"infer\"].",
        default="infer",
    )

    parser.add_argument(
        "--mean",
        type=float,
        help="For a given isoform, the mean tpm/count across samples \
        must be at least this to be considered, else discard isoform. [default: 0.0]",
        default=0.0,
    )

    parser.add_argument(
        "--median",
        type=float,
        help="For a given isoform, the median tpm/count across samples \
        must be at least this to be considered, else discard isoform. [default: 0.0]",
        default=0.0,
    )

    parser.add_argument(
        "--min",
        type=float,
        help="For a given isoform, each sample must have at least \
        this tpm/count to be considered, else discard isoform. [default: 0.0]",
        default=0.0,
    )

    parser.add_argument(
        "--max",
        type=float,
        help="For a given isoform, at least one sample must have at least \
        this tpm/count to be considered, else discard isoform. [default: 0.0]",
        default=0.0,
    )

    parser.add_argument(
        "-n", "--top_n",
        type=int,
        help="For a given gene, rank all isoforms by median_tpm, keep the \
        top_n ranked isoforms and discard the rest. '-1' to keep all isoforms. [default: -1]",
        default=-1,
    )

    parser.add_argument(
        "-c", "--top_n_column",
        type=str,
        choices=["median", "mean", "min", "max"],
        help="The metric by which to rank and filter top isoforms. \
            Used with '-n' to select expressed transcripts. [default: median]",
        default="median",
    )

    parser.add_argument(
        "-m", "--model",
        choices=["metagene", "consensus", "longest", "shortest"],
        help="Whether to output 'metagene', 'consensus', 'longest', 'shortest' model. \
        'longest' and 'shortest' select, for a given gene, the transcript with the longest \
        or shortest CDS, for now noncoding genes are ignored. Output is always after tpm \
        filtering has been applied. Multiple entries are allowed e.g. --model metagene \
        consensus longest [default: None]",
        type=str,
        default=[],
        nargs='*'
    )

    parser.add_argument(
        "-w ", "--tss_window",
        nargs=2,
        type=int,
        default=None,
        help="Pass two, space-separated, integers to specifiy the bp window around the TSS \
            as '<upstream>' '<downstream>'. Strand-orientation is inferred, i.e. '<upstream>' \
            will be in the 5\' direction of the TSS and <downstream> in the 3\' \
            direction. e.g. --tss_window 250 150. [default: None]"
    )

    parser.add_argument(
        "-e ", "--tes_window",
        nargs=2,
        type=int,
        default=None,
        help="Pass two, space-separated, integers to specifiy the bp window around the \
            transcription end site, TES, as '<upstream>' '<downstream>'. Strand-orientation \
            is inferred, i.e. '<upstream>' will be in the 5\' direction of the TES and \
            <downstream> in the 3\' direction. e.g. --tss_window 0 150. [default: None]"
    )

    parser.add_argument(
        "-x ", "--tx_to_gene",
        type=str,
        default=None,
        help="A TSV with transcript IDs in the first column and Gene IDs in the second. \
        The transcript IDs must match the first column entries of the --quant_files. \
        If this is not provided it will be deduced from the GTF/GFF3 and saved as \
        './annotations/intermediateFiles/tx2gene.tsv'."
    )

    args = parser.parse_args()
    print('\n')

        # List of arguments to check if they are provided without '--tpm_files'
    tpm_dependent_args = ['mean', 'median', 'min', 'max', 'top_n', 'top_n_column']

    # Check if any tpm_dependent_args are set without '--tpm_files'
    if not args.tpm_files:
        for arg in tpm_dependent_args:
            # Check if the argument is not its default value
            if getattr(args, arg) != parser.get_default(arg):
                print(f"\t\tWarning: Argument --{arg} is provided without --tpm_files. This will be ignored.")

    return args

def save_tss_tes_bed(args, GTF_path, GTF_file):
    """
    Generates and saves a BED file with Transcription Start Sites (TSS) from a GTF file.

    Parameters:
    - args: An object containing the command line arguments. Expected to have 'model'
      and 'tss_window' attributes.
    - GTF_path (str): The path to the input GTF file.
    - GTF_file (str): The filename of the GTF file.

    Notes:
    - The function assumes the existence of the './annotations/' directory.
    - The generated BED file contains TSS positions for each transcript in the GTF file.
    """

    # hacky approach to allow TSS or TES not be output
    output_TSS, output_TES = True, True
    if not args.tss_window:
        args.tss_window = [1, 1]
        output_TSS = False
    if not args.tes_window:
        args.tes_window = [1, 1]
        output_TES = False

    tss_bed, tes_bed = gtf_to_tss_tes_bed(GTF_path, tss_upstream=args.tss_window[0],
                                        tss_downstream=args.tss_window[1],
                                        tes_upstream=args.tes_window[0],
                                        tes_downstream=args.tes_window[1])

    tss_bed_out = './annotations/TSS_' + '.'.join(GTF_file.split('.')[:-1]) + '.bed'
    tes_bed_out = './annotations/TES_' + '.'.join(GTF_file.split('.')[:-1]) + '.bed'

    if output_TSS:
        print('\n\tSaving TSS:\t' + tss_bed_out + '\n')
        with open(tss_bed_out, 'w') as f:
            for entry in tss_bed:
                f.write(entry + '\n')

    if output_TES:
        print('\n\tSaving TES:\t' + tes_bed_out + '\n')
        with open(tes_bed_out, 'w') as f:
            for entry in tes_bed:
                f.write(entry + '\n')

def main():

    args = parse_arguments()

    create_output_directory(base_dir="./annotations/",output_prefix="intermediateFiles/")

    GTF_path = args.gtf
    GTF_file = args.gtf.split('/')[-1]

    if args.tpm_files:

        print('\tProcessing isoform quantification files')

        print('\n\tRemoving transcripts below threshold')

        transcript_df = filter_dataframe(process_files(args.tpm_files), \
                                        tscript_max = args.max, \
                                        tscript_min = args.min, \
                                        tscript_mean = args.mean, \
                                        tscript_median = args.median )

        tx_to_gene = {}

        if not args.tx_to_gene:
            print('\n\tGenerating transcript-gene relationships')
            tx_to_gene = extract_transcript_gene_relationship(GTF_path)
            tx2gene = './annotations/intermediateFiles/tx2gene.tsv'
            print('\n\tSaving transcript-gene relationships to:\t' + tx2gene)
            with open("./annotations/intermediateFiles/tx2gene.tsv", 'w') as f:
                for key, value in tx_to_gene.items():
                    f.write(f"{key}\t{value}\n")
        else:
            print('\n\tImporting transcript-gene relationships from ' + args.tx_to_gene)
            with open(args.tx_to_gene, 'r') as f:
                for line in f:
                    key, value = line.strip().split('\t')
                    tx_to_gene[key] = value

        result_df = add_gene_ids_and_subset(transcript_df, relationship=tx_to_gene, col = args.top_n_column, n = args.top_n)
        result_df_out = './annotations/intermediateFiles/filtered_' + '.'.join(GTF_file.split('.')[:-1]) + '.tsv'

        print('\tSaving quantification file to:\t\t' + result_df_out)
        result_df.to_csv(result_df_out, sep="\t", index=False)

        #these are the transcript ids that should be extracted from GTF/GFF
        transcript_ids = set(result_df['transcript_id'].tolist())

        #transcript_ids = extract_transcript_ids_from_dataframe(result_df)
        filtered_gtf_lines = filter_gtf_by_transcript_ids(input_file = GTF_path, transcript_ids = transcript_ids)

        GTF_file = 'filtered_' + '.'.join(GTF_file.split('.')[:-1]) + '.gtf'
        GTF_path = './annotations/' + GTF_file
        print('\tSaving transcript filtered GTF to:\t' + GTF_path)
        with open(GTF_path, 'w') as f:
            for line in filtered_gtf_lines:
                f.write(line)

        # Save TSS file for filtered GTF
        if args.tss_window or args.tes_window:
            save_tss_tes_bed(args, GTF_path = GTF_path, GTF_file = GTF_file)

    if args.model:
        cur_GTF_file = GTF_file
        cur_GTF_path = GTF_path

        if "shortest" or "longest" in args.model:
            longest_df, shortest_df = parse_gtf_for_cds_extremes(cur_GTF_path)

            if "longest" in args.model:
                GTF_file = 'longestCDS_' + '.'.join(cur_GTF_file.split('.')[:-1]) + '.gtf'
                GTF_path = './annotations/' + GTF_file
                print('\n\tSaving longest CDS GTF to: ' + GTF_path)

                #longest_df.drop(['gene_id', 'transcript_id', 'cds_length'], axis=1).to_csv(GTF_path, sep='\t', index=False, header = False, quoting=3)
                longest_df.to_csv(GTF_path, sep='\t', index=False, header = False, quoting=3)

                if args.tss_window or args.tes_window:
                    save_tss_tes_bed(args, GTF_path = GTF_path, GTF_file = GTF_file)

            if "shortest" in args.model:
                GTF_file = 'shortestCDS_' + '.'.join(cur_GTF_file.split('.')[:-1]) + '.gtf'
                GTF_path = './annotations/' + GTF_file
                print('\tSaving shortest CDS GTF to: ' + GTF_path)

                #shortest_df.drop(['gene_id', 'transcript_id', 'cds_length'], axis=1).to_csv(GTF_path, sep='\t', index=False, header=False, quoting=3)
                shortest_df.to_csv(GTF_path, sep='\t', index=False, header = False, quoting=3)
                if args.tss_window or args.tes_window:
                    save_tss_tes_bed(args, GTF_path = GTF_path, GTF_file = GTF_file)


        if "metagene" in args.model:
            output_str = create_metagene_model(cur_GTF_path)

            GTF_file = 'meta_' + '.'.join(cur_GTF_file.split('.')[:-1]) + '.gtf'
            GTF_path = './annotations/' + GTF_file
            print('\tSaving metagene GTF to: ' + GTF_path)
            with open(GTF_path, 'w') as f:
                f.write(output_str)

            if args.tss_window or args.tes_window:
                save_tss_tes_bed(args, GTF_path = GTF_path, GTF_file = GTF_file)

        if "consensus" in args.model:

            output_str, genes_without_consensus = create_constitutive_model(cur_GTF_path)

            GTF_file = 'consensus_' + '.'.join(cur_GTF_file.split('.')[:-1]) + '.gtf'
            GTF_path = './annotations/' + GTF_file

            print('\tSaving consensus GTF to: ' + GTF_path)
            with open(GTF_path, 'w') as f:
                f.write(output_str)

            genes_without_consensus_out = './annotations/genes_without_consensus_model.txt'
            print('\n\tA CONSENSUS MODEL COULD NOT BE GENERATED FOR ' + str(len(genes_without_consensus)) + ' GENES')
            print('\tIf this number is large, consider filtering by TPM expression more strictly or using a more conservative GTF.')
            print('\tIf this number is small, consider manually removing problematic transcripts from the quantification TSVs and rerunning this module.')
            print('\tSaving genes for which there is no consensus model to:\t' + genes_without_consensus_out)
            with open(genes_without_consensus_out, 'w') as f:
                f.writelines(gene_id + "\n" for gene_id in genes_without_consensus)

            if args.tss_window or args.tes_window:
                save_tss_tes_bed(args, GTF_path = GTF_path, GTF_file = GTF_file)


    elif args.model != "none" and ( args.tss_window or args.tes_window):
        save_tss_tes_bed(args, GTF_path = GTF_path, GTF_file = GTF_file)

    print('\n')


if __name__ == "__main__":
    main()
