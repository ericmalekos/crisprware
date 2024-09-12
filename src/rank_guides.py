#!/usr/bin/env python3

'''
    This script takes the output from score_guides and ranks them according to
    user defined criteria
'''

import pandas as pd
import argparse
import pybedtools as pyb
import matplotlib.pyplot as plt
#from utils.utility_functions import create_output_directory
from utils.rank_guides_functions import create_combined_weighted_column,\
    validate_and_modify_bed,df_to_pybed,gRNA_to_bed,select_guides,\
        gRNA_to_tscript,group_and_minimize,analyze_target_ids
from utils.gtf_bed_processing_functions import truncate_gtf
from utils.utility_functions import create_output


def restricted_int(x):
    x = int(x)
    if x < 0 or x > 100:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return x

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Rank gRNA output from score_guides."
    )

    parser.add_argument(
        "-k", "--scored_guides",
        type=str,
        help="<score_guides_output>.bed output from score_guides.",
        required=True
    )

    parser.add_argument(
        "-t", "--targets",
        type=str,
        help="BED/GTF/GFF used to select final guides per target. \
            For GTF/GFF, set --target_mode to either 'gene' or 'transcript'. \
            For BED, targets are each entry. Use '--number_of_targets' \
            to set the number of guides chosen for each target.",
        required=True
    )

    parser.add_argument(
        "--target_mode",
        type=str,
        default="gene",
        help="If a GTF/GFF is used to select targets, gRNAs can be grouped \
        at either the 'tx' or 'gene' level e.g. '--target_mode gene -n 10' \
        chooses 10 guides per gene, '--target_mode tx -n 10' chooses 10 \
        per transcript [default: gene].",
        choices=["gene", "tx"]
    )

    parser.add_argument(
        "-f", "--feature",
        type=str,
        help="If GTF/GFF passed, use this feature for processing e.g. 'exon', 'CDS', '5UTR', etc. \
            The feature appears in the third column of the GTF/GFF [default: CDS]. ",
        default='CDS'
    )

    parser.add_argument(
        "-p", "--percentile_range",
        nargs=2,
        help="Allowable range of guide for each transcript and feature set, e.g. \n \
            '-p 60 80 -f exon' returns gRNAs in the 60th to 80th percentile of exons \n \
            for a given transcript. Default setting returns guides \n \
            anywhere in the CDS for each transcript [default: 0 100]",
        type=restricted_int,
        default=[0, 100]
    )

    parser.add_argument(
        "-n", "--number_of_guides",
        type=int,
        help="Number of guides returned per target.'-1' to keep all guides [default: -1]",
        default=-1
    )
    parser.add_argument(
        "--min_spacing",
        type=int,
        help="The minimum nucleotide space between guides for a given target. \
            e.g. --min_spacing 10, requires guides 10 nts appart. 0 to allow overlapping guides.[default: 0]",
        default=0
    )

    parser.add_argument(
        "--output_all",
        default=False,
        action="store_true",
        help="Set flag to save gRNA-target TSVs at each \
            stage of filtering rather than just the end.[default: False]"
    )

    parser.add_argument(
        "--plot_histogram",
        default=False,
        action="store_true",
        help="Set flag to plot a histogram of the distribution of gRNAs per target \
            after each filtering step. Sets '--output_all' to True.[default: False]"
    )

    parser.add_argument(
        "-c", "--filtering_columns",
        type=str,
        help="One or more space-separated column names used for filtering. Uses raw values. \
            e.g. '-c rs3_z_score specificity_Hg38_index'.",
        required=False,
        nargs='*'
    )

    parser.add_argument(
        "-m", "--minimum_values",
        type=float,
        help="A space-separated list of minimum values for each column in passed \
        by --ranking_columns. e.g. '-c rs3_z_score specificity_Hg38_index -m \"-1\" 0.2' \
        Default is no minimum [default: None]",
        required=False,
        nargs='*'
    )

    parser.add_argument(
        "-r", "--ranking_columns",
        type=str,
        help="One or more space-separated column names used for guide ranking. \
            e.g. '-r rs3_score_Hsu2013 rs3_score_Chen2013'.",
        required=False,
        nargs='*'
    )

    parser.add_argument(
        "-w", "--column_weights",
        type=float,
        help="A space-separated list of weight values for each column in passed \
        by --ranking_columns. e.g. '-c rs3_score specificity_Hg38_index -w 1 0' \
        Default is equal weighting for all ranking columns.",
        required=False,
        nargs='*'
    )

    parser.add_argument(
        "--normalize_columns",
        help="Scale ranking column values to 0 to 1 [default: True]",
        default=True,
        action="store_false"
    )

    parser.add_argument(
        "-o", "--output_directory",
        help="Path to output. [default: current directory]",
        type=str,
        default=""
    )

    return parser.parse_args()


def validate_scored_gRNA_columns(header, columns):
    required_columns = ['#chr', 'start', 'stop', 'strand']
    # Check the fixed columns
    for column in columns:
        if column not in header:
            raise ValueError("\t\"" + column + "\" column not present")

    for column in required_columns:
        if column not in header:
            raise ValueError("\t" + column + "not present, check input")



def update_no_targets(original_target_ids, target_df):
    updated_target_ids = set(target_df['target_id'].unique())
    removed_target_ids = original_target_ids - updated_target_ids
    return removed_target_ids

def create_histogram(df, output_pdf_path, target_type="Target ID", targetsWith0Guides=""):
    """
    Creates and saves a histogram that shows the frequency of gRNAs per target ID,
    including targets with zero gRNAs, to a PDF file.

    Parameters:
    - df (pd.DataFrame): DataFrame containing gRNA data, including a 'target_id' column.
    - output_pdf_path (str): Path (including filename) where the histogram PDF will be saved.
    - targetsWith0Guides (list): A list of target IDs that have zero associated gRNAs.

    The histogram includes an additional bar representing the number of targets with zero gRNAs.
    """
    # Count the occurrences of each 'target_id'
    target_id_counts = df['target_id'].value_counts()

    # Count the frequencies of these occurrences
    frequency_counts = target_id_counts.value_counts().sort_index()

    # If there are targets with 0 guides, add this count to the frequency_counts
    if targetsWith0Guides:
        frequency_counts[0] = len(targetsWith0Guides)

    frequency_counts = frequency_counts.sort_index()

    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    ax = frequency_counts.plot(kind='bar')
    plt.title(f'Frequency of {target_type} counts')
    plt.xlabel(f'Number of gRNAs per {target_type}')
    plt.ylabel(f'Number of unique {target_type}s')

    # Setting x-axis labels to only appear every 5 ticks
    ticks = ax.xaxis.get_ticklocs()
    tick_labels = [str(int(label)) for label in frequency_counts.index]  # Convert index to string for tick labels
    ax.xaxis.set_ticks(ticks[::5])
    ax.xaxis.set_ticklabels(tick_labels[::5])

    plt.tight_layout()

    # Save the figure to a PDF file
    plt.savefig(output_pdf_path)
    plt.close()

def save_notarget_set(noSgRNATargetSet, filename):

    sorted_strings = sorted(noSgRNATargetSet)
    with open(filename, "w+") as file:
        for item in sorted_strings:
            file.write(item + "\n")


def filter_df_by_cutoff(df, column_name, cutoff_value):
    """
    Parameters:
    df (pd.DataFrame): The DataFrame to filter.
    column_names (list): A column name
    cutoff_values (list): A minimal column value

    Returns:
    pd.DataFrame, dict: The filtered DataFrame and a dictionary with column names as keys and lists of removed target_ids as values.
    """
    df = df[df[column_name] >= cutoff_value]

    print(f"\tRemaining gRNAs after applying cutoff {cutoff_value} to '{column_name}': {len(df)}")

    return df

def save_output(args,finalgRNAs,out_path,initial_target_ids):

    targetsWithoutgRNAs = update_no_targets(initial_target_ids,finalgRNAs)

    _ = analyze_target_ids(finalgRNAs, targetsWithoutgRNAs)

    if args.output_all: finalgRNAs.to_csv(out_path, sep='\t', index=False)
    if args.plot_histogram: create_histogram(finalgRNAs, out_path[:-3] + "pdf", args.target_mode, targetsWithoutgRNAs)

def main():

    pd.set_option('display.max_columns', None)
    args = parse_arguments()

    if args.ranking_columns and not args.column_weights:
        print(f"\n\n\tColumn weights not set, setting to weights to 1")
        args.column_weights = [1] * len(args.ranking_columns)

    if args.ranking_columns and (len(args.ranking_columns) != len(args.column_weights)):
        raise ValueError("\n\n\tThe length of --ranking_columns and --minimum_values must match.\n\n")

    if args.filtering_columns and (not args.minimum_values or (len(args.filtering_columns) != len(args.minimum_values))):
        raise ValueError("\n\n\tThe length of --filtering_columns and --minimum_values must match.\n\n")

    gRNA_output_path, tmp_path = create_output(args.scored_guides, outdir=args.output_directory, extension="rankedgRNA", stripped="_scoredgRNA", tmp=True)
    gRNA_output_path += ".bed"

    with open(args.scored_guides, 'r') as file: header = file.readline().strip().split('\t')


    if args.ranking_columns:
        validate_scored_gRNA_columns(header,columns = args.ranking_columns)
    if args.filtering_columns:
        validate_scored_gRNA_columns(header,columns = args.filtering_columns)

    gRNAs = pyb.BedTool(args.scored_guides)

    # Determine target file type
    targetFileType = args.targets.split(".")[-1].lower()
    initial_target_ids = None
    targetsWithoutgRNAs = None


    if targetFileType in ["gff", "gtf", "gff3", "gff2"]:
        print(f'\n\t{args.targets} is {targetFileType.upper()} format')

        # For establishing initial count
        _, transcript_ids, gene_ids = truncate_gtf(input_file = args.targets, feature = args.feature, percentiles = [0,100])

        finalgRNAs, targetIDs = gRNA_to_tscript(gRNAs=gRNAs, mode=args.target_mode, targets=pyb.BedTool(args.targets), header=header)
        initial_target_ids = transcript_ids if args.target_mode == "tx" else gene_ids
        finalgRNAs = finalgRNAs.drop_duplicates()
        print(f'\n\tInitial {args.target_mode} count:\t{len(initial_target_ids)}')

        print("\n\n\tPrior to positional filtering: ")
        out_path = tmp_path + "0_gRNAsTargets_preFilter.bed"
        save_output(args,finalgRNAs,out_path,initial_target_ids)

        trimmed = tmp_path + "trimmedAnnotation_" + args.targets.split("/")[-1]

        trimmed_gtf, transcript_ids, gene_ids = truncate_gtf(input_file = args.targets, feature = args.feature, percentiles = args.percentile_range)

        if args.target_mode == "tx": initial_target_ids = transcript_ids
        else: initial_target_ids = gene_ids
        #print(f'\n\tTarget count after {args.feature} trimmming:\t{len(initial_target_ids)}');

        trimmed_gtf.to_csv(trimmed, sep="\t", header=False, index=False, doublequote=False, quoting=3, quotechar="",  escapechar="\\")
        targets = pyb.BedTool(trimmed)
        finalgRNAs, targetIDs = gRNA_to_tscript(gRNAs=gRNAs, mode=args.target_mode, targets=targets, header=header)

    else:
        if targetFileType not in ["bed"]:
            print(f'{args.targets} file extension unrecognized, for expected behavior pass GFF/GTF/BED file.\n \
            Defaulting to BED behavior')

        target_bed_df, initial_target_ids = validate_and_modify_bed(args.targets)
        print(f"\n\tNumber of targets: {len(initial_target_ids)}")
        target_bed, target_header = df_to_pybed(df=target_bed_df)

        finalgRNAs,targetIDs = gRNA_to_bed(gRNAs=gRNAs, targets=target_bed, header=header, target_header=target_header)


    if args.ranking_columns:
        finalgRNAs = create_combined_weighted_column(df=finalgRNAs, column_names=args.ranking_columns,
                                                  weights=args.column_weights, normalize_columns=args.normalize_columns)


    targetsWithoutgRNAs = initial_target_ids - targetIDs

    # filter by column values
    print("\n\n\tPrior to filtering: ")
    out_path = tmp_path + "1_gRNAsTargets_preFilter.bed"
    save_output(args=args,finalgRNAs=finalgRNAs,out_path=out_path,initial_target_ids=initial_target_ids)

    i = 2
    if args.filtering_columns:
        for column, cutoff in zip(args.filtering_columns, args.minimum_values):
            print(f"\n\tFiltering column {column} > {cutoff}")
            finalgRNAs = filter_df_by_cutoff(df = finalgRNAs, column_name=column, cutoff_value=cutoff)
            out_path = tmp_path + str(i) + "_gRNAsTargets_" + column + str(cutoff) + "Filter.bed"
            save_output(args=args,finalgRNAs=finalgRNAs,out_path=out_path,initial_target_ids=initial_target_ids)
            i += 1

    # Filter by spacing
    out_path = tmp_path + str(i) + "_gRNAsTargets_spacing" + "Filter.bed" ; i += 1
    if args.ranking_columns:
        finalgRNAs = select_guides(df=finalgRNAs, rank_column="combined_weighted", min_spacing=args.min_spacing)
        save_output(args,finalgRNAs,out_path,initial_target_ids)

    # Filter by number_of_guides per target
    #if args.number_of_guides != -1:
    if args.ranking_columns:
        finalgRNAs = group_and_minimize(finalgRNAs, "combined_weighted", args.number_of_guides)

    counts = analyze_target_ids(finalgRNAs, update_no_targets(initial_target_ids,finalgRNAs))

    if 'id,sequence,pam,chromosome,position,sense' in finalgRNAs.columns:
        finalgRNAs[['sequence']] = finalgRNAs['id,sequence,pam,chromosome,position,sense'].str.split(',', expand=True).iloc[:, 1:2]
        finalgRNAs = finalgRNAs.drop(['id,sequence,pam,chromosome,position,sense'], axis = 1)

    if 'context' in finalgRNAs.columns:
        finalgRNAs = finalgRNAs.drop(['context'], axis = 1)

    finalgRNAs.to_csv(gRNA_output_path, sep='\t', index=False)
    if args.plot_histogram: create_histogram(finalgRNAs, gRNA_output_path[:-3] + "pdf", args.target_mode, update_no_targets(initial_target_ids,finalgRNAs))

    save_notarget_set(noSgRNATargetSet = targetsWithoutgRNAs.union(update_no_targets(initial_target_ids,finalgRNAs)),
                      filename=gRNA_output_path.replace("rankedgRNA.bed", "") + "targetsWith0gRNAs.txt")

    counts = counts.reset_index()
    counts.columns = ['target_id', 'gRNA_count']
    counts.to_csv(gRNA_output_path.strip("rankedgRNA.bed") + "targetgRNACounts.bed", sep='\t', index=False)


if __name__ == "__main__":
    main()
