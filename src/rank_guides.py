#!/usr/bin/env python3

'''
    This script takes the output from score_guides and ranks them according to
    user defined criteria 
'''

import pandas as pd
#from pathlib import Path
import argparse
import pybedtools as pyb
import matplotlib.pyplot as plt
from utils.utility_functions import create_output_directory
from utils.rank_guides_functions import create_combined_weighted_column,\
    validate_and_modify_bed,df_to_pybed,sgRNA_to_bed,select_guides,\
        sgRNA_to_tscript,group_and_minimize,analyze_target_ids
from utils.gtf_bed_processing_functions import truncate_gtf

#


def restricted_int(x):
    x = int(x)
    if x < 0.0 or x > 100.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0, 100]")
    return x

def restricted_float(x):
    x = int(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return x

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Rank sgRNA output from score_guides."
    )

    parser.add_argument(
        "-k", "--scored_guides",
        type=str,
        help="<score_guides_output>.tsv output from score_guides.",
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
        help="If a GTF/GFF is used to select targets, sgRNAs can be grouped \
        at either the 'transcript' or 'gene' level e.g. '--target_mode gene -n 10' \
        chooses 10 guides per gene, '--target_mode transcript -n 10' chooses 10 \
        per transcript [default: gene].",
        choices=["gene", "transcript"]
    )

    parser.add_argument(
        "-f", "--feature",
        type=str,
        help="If GTF/GFF passed, use this feature for processing e.g. 'exon', 'CDS', '5UTR', etc. \
            The feature appears in the third column of the GTF/GFF.",
        default='CDS'
    )

    parser.add_argument(
        "-p", "--percentile_range",
        help="Allowable range of guide for each transcript and feature set, e.g. \n \
            '-p 60 -f exon' returns sgRNAs in the 5'-most 60th percentile of exons \n \
            for a given transcript. Default setting returns guides \n \
            anywhere in the CDS for each transcript [default: 100]",
        type=restricted_int,
        default=100
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
        help="The minimum nucleotide space between guides for a single target. \
            0 to allow overlapping guides.[default: 10]", 
        default=10
    )

    parser.add_argument(
        "--output_all", 
        default=False,
        action="store_true", 
        help="Set flag to save sgRNA-target TSVs at each \
            stage of filtering rather than just the end.[default: False]"
    )
    
    parser.add_argument(
        "--plot_histogram", 
        default=False,
        action="store_true", 
        help="Set flag to plot a histogram of the distribution of sgRNAs per target \
            after each filtering step. Sets '--output_all' to True.[default: False]"
    )    

    parser.add_argument(
        "-c", "--ranking_columns",
        type=str,
        help="One or more space-separated column names used for guide ranking. \
            e.g. '-c rs3_z_score specificity_Hg38_index'.",
        required=True,
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
        "-w", "--column_weights",
        type=float,
        help="A space-separated list of weight values for each column in passed \
        by --ranking_columns. e.g. '-c rs3_z_score specificity_Hg38_index -m \"-1\" 0.2' \
        Default is 1 for all columns, i.e. equal weighting [default: 1]",
        required=False,
        nargs='*'
    )

    parser.add_argument(
        "-o", "--output_prefix", 
        type=str, 
        help="Prefix for output file", 
        default=""
    )

    return parser.parse_args()

    
def validate_scored_sgRNA_columns(header, rank_columns):
    required_columns = ['#chr', 'start', 'stop', 'sequence', 'strand']
    # Check the fixed columns
    for column in rank_columns: 
        if column not in header:
            raise ValueError("\t\"" + column + "\" column not present")
        
    for column in required_columns:
        if column not in header:
            raise ValueError("\t" + column + "not present, check input")


def update_no_targets(original_target_ids, target_df):
    updated_target_ids = set(target_df['target_id'].unique())
    removed_target_ids = original_target_ids - updated_target_ids
    return removed_target_ids

def create_histogram(df, output_pdf_path, targetsWith0Guides=""):
    """
    Creates and saves a histogram that shows the frequency of sgRNAs per target ID, 
    including targets with zero sgRNAs, to a PDF file.

    Parameters:
    - df (pd.DataFrame): DataFrame containing sgRNA data, including a 'target_id' column.
    - output_pdf_path (str): Path (including filename) where the histogram PDF will be saved.
    - targetsWith0Guides (list): A list of target IDs that have zero associated sgRNAs.
    
    The histogram includes an additional bar representing the number of targets with zero sgRNAs.
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
    plt.title('Frequency of Target ID Counts')
    plt.xlabel('Number of sgRNAs per Target ID')
    plt.ylabel('Number of Unique Target IDs')

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
    with open(filename, "w") as file:
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
        
    print(f"\tRemaining sgRNAs after applying cutoff {cutoff_value} to '{column_name}': {len(df)}")
    
    return df

def save_output(args,finalsgRNAs,out_path,initial_target_ids):
    
    targetsWithoutsgRNAs = update_no_targets(initial_target_ids,finalsgRNAs)

    _ = analyze_target_ids(finalsgRNAs, targetsWithoutsgRNAs)

    if args.output_all: finalsgRNAs.to_csv(out_path, sep='\t', index=False)
    if args.plot_histogram: create_histogram(finalsgRNAs, out_path[:-3] + "pdf", targetsWithoutsgRNAs)

def main():
    pd.set_option('display.max_columns', None)
    args = parse_arguments()

    if len(args.ranking_columns) != len(args.minimum_values):
        raise ValueError("The length of --ranking_columns and --minimum_values must match.")

    sgRNA_output_path = "./" + args.output_prefix + "RankedSgRNAs/" + args.output_prefix.split("/")[-1] + "RankedSgRNA.tsv"
    tmp_path = create_output_directory(base_dir="./" + args.output_prefix + "RankedSgRNAs/",output_prefix="tmp/")

    with open(args.scored_guides, 'r') as file: header = file.readline().strip().split('\t')        
    validate_scored_sgRNA_columns(header,rank_columns = args.ranking_columns)

    sgRNAs = pyb.BedTool(args.scored_guides)
    
    # Determine target file type
    targetFileType = args.targets.split(".")[-1].lower()
    initial_target_ids = None
    targetsWithoutsgRNAs = None


    if targetFileType in ["gff", "gtf", "gff3", "gff2"]:
        print(f'\n\t{args.targets} is {targetFileType.upper()} format')


        _, transcript_ids, gene_ids = truncate_gtf(input_file = args.targets, feature = args.feature, percentile = 100)

        finalsgRNAs, targetIDs = sgRNA_to_tscript(sgRNAs=sgRNAs, mode=args.target_mode, targets=pyb.BedTool(args.targets), header=header)
        initial_target_ids = transcript_ids if args.target_mode == "transcript" else gene_ids
        print(f'\n\tInitial {args.target_mode} count:\t{len(initial_target_ids)}')

        print("\n\n\tPrior to filtering: ")
        out_path = tmp_path + args.output_prefix + "0_sgRNAsTargets_preFilter.tsv"
        save_output(args,finalsgRNAs,out_path,initial_target_ids)

        trimmed = tmp_path + "trimmedAnnotation_" + args.targets.split("/")[-1]

        trimmed_gtf, transcript_ids, gene_ids = truncate_gtf(input_file = args.targets, feature = args.feature, percentile = args.percentile_range)
        
        if args.target_mode == "transcript": initial_target_ids = transcript_ids
        else: initial_target_ids = gene_ids
        print(f'\n\tTarget count after CDS trimmming:\t{len(initial_target_ids)}'); 

        trimmed_gtf.to_csv(trimmed, sep="\t", header=False, index=False, doublequote=False, quoting=3, quotechar="",  escapechar="\\")
        targets = pyb.BedTool(trimmed)
        finalsgRNAs, targetIDs = sgRNA_to_tscript(sgRNAs=sgRNAs, mode=args.target_mode, targets=targets, header=header)

    else:
        if targetFileType not in ["bed"]:
            print(f'{args.targets} file extension unrecognized, for expected behavior pass GFF/GTF/BED file.\n \
            Defaulting to BED behavior')
        
        target_bed_df, initial_target_ids = validate_and_modify_bed(args.targets)
        print(f"\n\tNumber of targets: {len(initial_target_ids)}")  
        target_bed, target_header = df_to_pybed(df=target_bed_df)

        finalsgRNAs,targetIDs = sgRNA_to_bed(sgRNAs=sgRNAs, targets=target_bed, header=header, target_header=target_header)

    finalsgRNAs = create_combined_weighted_column(df=finalsgRNAs, column_names=args.ranking_columns, weights=args.column_weights)

    targetsWithoutsgRNAs = initial_target_ids - targetIDs
    
    # filter by column values
    print("\n\n\tPrior to filtering: ")
    out_path = tmp_path + args.output_prefix + "1_sgRNAsTargets_preFilter.tsv"
    save_output(args=args,finalsgRNAs=finalsgRNAs,out_path=out_path,initial_target_ids=initial_target_ids)

    i = 2
    for column, cutoff in zip(args.ranking_columns, args.minimum_values):
        finalsgRNAs = filter_df_by_cutoff(df = finalsgRNAs, column_name=column, cutoff_value=cutoff)
        out_path = tmp_path + args.output_prefix + str(i) + "_sgRNAsTargets_" + column + str(cutoff) + "Filter.tsv"
        save_output(args=args,finalsgRNAs=finalsgRNAs,out_path=out_path,initial_target_ids=initial_target_ids)
        i += 1
    
    # Filter by spacing
    out_path = tmp_path + args.output_prefix + str(i) + "_sgRNAsTargets_spacing" + "Filter.tsv" ; i += 1
    finalsgRNAs = select_guides(df=finalsgRNAs, rank_column="combined_weighted", min_spacing=args.min_spacing)
    save_output(args,finalsgRNAs,out_path,initial_target_ids)

    # Filter by number_of_guides per target
    finalsgRNAs = group_and_minimize(finalsgRNAs, "combined_weighted", args.number_of_guides)
    counts = analyze_target_ids(finalsgRNAs, update_no_targets(initial_target_ids,finalsgRNAs))

    finalsgRNAs.to_csv(sgRNA_output_path, sep='\t', index=False)
    if args.plot_histogram: create_histogram(finalsgRNAs, sgRNA_output_path[:-3] + "pdf", update_no_targets(initial_target_ids,finalsgRNAs))

    save_notarget_set(noSgRNATargetSet = targetsWithoutsgRNAs.union(update_no_targets(initial_target_ids,finalsgRNAs)),
                      filename="."+sgRNA_output_path.strip("RankedSgRNA.tsv") + "TargetsWith0sgRNAs.tsv")
    
    counts.to_csv("." + sgRNA_output_path.strip("RankedSgRNA.tsv") + "TargetsgRNACounts.tsv", sep='\t', index=True)


if __name__ == "__main__":
    main()