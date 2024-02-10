#!/usr/bin/env python3

'''
    This script takes the output from FiltersgRNABED and scores it
    with Ruleset3 (cleavage efficieny, sequence-based) and Guidescan2 (off-target). 
'''

import pandas as pd
from pathlib import Path
import argparse
import pybedtools as pyb
import re
import matplotlib.pyplot as plt


def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return x

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Rank sgRNAs from ScoresgRNAs."
    )

    parser.add_argument(
        "-k", "--scoredGuides",
        type=str,
        help="AllScoredGuides.tsv output from ScoresgRNAs.",
        required=True
    )

    parser.add_argument(
        "-t", "--targets",
        type=str,
        help="BED/GTF/GFF used to select final guides per target. \
            For GTF/GFF, targets are transcripts. For BED, targets are each entry. \
            Use '--number_of_targets' to set the number of guides chosed for each target.",
        required=True
    )

    parser.add_argument(
        "--target_mode",
        type=str,
        default="gene",
        help="If a GTF/GFF is used to select targets, sgRNAs can be grouped \
        at either the 'transcript' or 'gene' level [default: gene].",
        choices=["gene", "transcript"]
    )

    parser.add_argument(
        "-f", "--feature",
        type=str,
        help="If GTF/GFF passed, use this feature for processing e.g. 'exon', 'CDS', '5UTR', etc. \
            The feature appear in the third column of the GTF/GFF.",
        default='CDS'
    )

    parser.add_argument(
        "-p", "--percentile_range",
        help="Allowable range of guide for each transcript and feature set, e.g. \n \
            '0.5' returns sgRNAs in the 5'-most half of the relevant \n \
            features for a transcript. Default settings returns guides \n \
            anywhere in the exons for each transcript [default: 1]",
        type=restricted_float,
        default=1
    )

    parser.add_argument(
        "-n", "--number_of_guides",
        type=int,
        help="Number of guides returned per target.'-1' to keep all guides [default: -1]",
        default=-1
    )

    parser.add_argument(
        "--min_specificity", 
        type=restricted_float, 
        default=0.2, 
        help="Minimum guidescan specificity score from [0,1]")

    parser.add_argument(
        "--min_cleavage", 
        type=float, 
        default=-1, 
        help="Minimum cleavage score after scaling \
            scores to [0,1]. Used with --score_column")

    parser.add_argument(
        "--specificity_column", 
        type=str, 
        help="Off-target specificty column to use with --min_specificity. \
            If not provided, will use the first column \
            with 'specificity' in the column name.", 
        default=None)

    parser.add_argument(
        "--cleavage_column", 
        type=str, 
        default="rs3_z_score", 
        help="Minimum cleavage score used with --min_cleavage. \
        Uses rs3_z_score as default but can be set to \
        any column that contains numeric entries.")
    
    parser.add_argument(
        "--rank_column", 
        type=str, 
        default="rs3_z_score", 
        help="Column used for final ranking and guide selection. \
            Set to 'combined_weighted' for the weighted average of \
            the --specificity_column and --cleavage_column.[default: 'rs3_z_score']")

    parser.add_argument(
        "--specificity_weight", 
        type=float, 
        default=1.0, 
        help="Weight for cleavage scoring")

    parser.add_argument(
        "--cleavage_weight", 
        type=float, 
        default=1.0, 
        help="Weight for guidescan")

    parser.add_argument(
        "--min_spacing", 
        type=int,
        help="The minimum nucleotide space between guides for a single target. \
            0 to allow overlapping guides.[default: 10]", 
        default=10)

    parser.add_argument(
        "--output_all", 
        default=False,
        action="store_true", 
        help="Set flag to save sgRNA-target TSVs at each \
            stage of filtering rather than just the end.[default: False]")
    
    parser.add_argument(
        "--plot_histogram", 
        default=False,
        action="store_true", 
        help="Set flag to plot a histogram of the distribution of sgRNAs per target \
            after each filtering step. Sets '--output_all' to True.[default: False]")    

    parser.add_argument(
        "--prefix", 
        type=str, 
        help="Prefix for output file", 
        default="")

    return parser.parse_args()


def extract_transcript_id(attributes):
    """
    Extracts the transcript_id from the provided attributes string. The function is designed to handle both 
    GTF and GFF file formats. In case of GTF format, the transcript_id is enclosed in double quotes and followed by a semicolon.
    In case of GFF format, the transcript_id is followed by an equals sign and doesn't have the double quotes and the semicolon.

    Parameters:
    attributes (str): A string containing the attributes of a feature in a GTF/GFF file.

    Returns:
    str: The transcript_id extracted from the attributes. If no transcript_id is found, the function will print an error 
    message and terminate the program.

    Example:
    >>> extract_transcript_id('transcript_id "ENST00000456328";')
    'ENST00000456328'
    >>> extract_transcript_id('transcript_id=ENST00000456328.2;')
    'ENST00000456328.2'
    """

    match = re.search(r'transcript_id[= ]"?([^";]*)"?', attributes)
    #match = re.search('transcript_id[ =]"?(.*?)"?;?', attributes)
    #match = re.search(r'gene_id[= ]"?([^";]*)"?', attributes)
    if match:
        return match.group(1)
    else:
        print("ERROR: --feature '"  + feature + "' missing 'transcript_id' attribute")
        exit(keep_kernel=True)
        
def extract_ids(attributes):
    
    #print(attributes)
    
    transcript_match = re.search(r'transcript_id[= ]"?([^";]*)"?', attributes)

    gene_match = re.search(r'gene_id[= ]"?([^";]*)"?', attributes)

    #print(transcript_match.group(1))
    
    transcript_id = transcript_match.group(1) if transcript_match else None
    gene_id = gene_match.group(1) if gene_match else None
    
    #print(transcript_id, gene_id)
    return transcript_id, gene_id

def process_gtf(input_file, feature = "exon", percentiles = [0,100]):
    """
    Processes a GTF/GFF file to extract and truncate specified features based on percentiles.

    Parameters:
    input_file (str): Path to the input GTF/GFF file.
    feature (str): The feature of interest. Default is 'exon'.
    percentiles (list): A list specifying the lower and upper percentile thresholds to truncate the features. Default is (0,50).

    Returns:
    DataFrame: A pandas DataFrame containing the selected and truncated features. 

    Raises:
    ValueError: If the specified feature is not found in the GTF/GFF file.

    The function reads the GTF/GFF file and selects the rows of the specified feature. 
    It then calculates the length of each feature and the cumulative length for each transcript.
    The entries where 'cumulative_percentile' exceeds 'max_percentile' are identified and 
    the 'end' or 'start' position of these entries are adjusted accordingly. 
    Finally, the function returns a DataFrame containing the selected and truncated features.

    Example:
    >>> df = process_gtf('test.gtf', 'CDS', (0, 50))
    """


    percentiles.sort()

    print(f'\tProcessing: \t{input_file} \n\tFeature: \t{feature} \n\tPercentile range: \t{percentiles}')
    
    df = pd.read_csv(input_file, sep='\t', header=None, comment='#')
    df.columns = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]

    # Select the rows for the specified feature
    df = df[df["feature"] == feature]
    
    if df.empty:
        print("ERROR: --feature '"  + feature + "' not found in GTF/GFF column 3")
        exit(keep_kernel=True)
        #        raise ValueError("ERROR: --feature '"  + feature + "' not found in GTF/GFF column 3")


    df_positive = df[df["strand"] == "+"].sort_values("start")
    df_negative = df[df["strand"] == "-"].sort_values("start", ascending=False)

    # Combine the sorted dataframes
    df = pd.concat([df_positive, df_negative])

    # Extract the transcript_id from the attributes column

    #print(df["attributes"].apply(extract_ids))
    #df["transcript_id"], df["gene_id"] = df["attributes"].apply(extract_ids)
    
    df[["transcript_id", "gene_id"]] = df["attributes"].apply(extract_ids).apply(pd.Series)


    #df["transcript_id"] = df["attributes"].apply(extract_transcript_id)
    
    print(f'\n\tNumber of {feature} entries before processing:{df.shape[0]}')
    # TODO check this
    # print(f'\n\tNumber of unique transcripts before processing:{sorted(df["transcript_id"].drop_duplicates())}')

        
    # Calculate the length of each feature
    df["length"] = abs(df["end"] - df["start"]) + 1

    # Group by the "transcript_id" and calculate the cumulative length and total length
    df["cumulative_length"] = df.groupby("transcript_id")["length"].cumsum()
    df["total_length"] = df.groupby("transcript_id")["length"].transform(sum)

    # Calculate the percentile of the cumulative length
    df["cumulative_percentile"] = df["cumulative_length"] / df["total_length"] * 100

    # Identify the entries where 'cumulative_percentile' exceeds 'max_percentile'
    overlaps = df[df["cumulative_percentile"] > percentiles[1]].groupby('transcript_id').first().reset_index()

    # Adjust the 'end' or 'start' position of overlapping entries
    for i, row in overlaps.iterrows():
        excess_percentile = row['cumulative_percentile'] - percentiles[1]
        excess_length = int(row['total_length'] * excess_percentile / 100)
        if row['strand'] == '+':
            df.loc[(df['transcript_id'] == row['transcript_id']) & (df['end'] == row['end']), 'end'] -= excess_length
        else:
            df.loc[(df['transcript_id'] == row['transcript_id']) & (df['start'] == row['start']), 'start'] += excess_length
        #print(df.to_string())
        df.loc[(df['transcript_id'] == row['transcript_id']), 'length'] = abs(df['end'] - df['start']) + 1
    # Recalculate the cumulative length and the cumulative percentile
    df["cumulative_length"] = df.groupby("transcript_id")["length"].cumsum()

    df["cumulative_percentile"] = (df["cumulative_length"] / df["total_length"] * 100).astype(int)
    
    
    df = df[df["cumulative_percentile"] <= percentiles[1]]
    
    print(f'\n\tNumber of {feature} entries after processing:\t{df.shape[0]}')
    # TODO check this
    #print(f'\n\tNumber of unique transcripts after processing:{sorted(df["transcript_id"].drop_duplicates())}')
    
    # Drop the extra columns
    df = df.drop(columns=["length", "cumulative_length", "total_length", "cumulative_percentile", "transcript_id", "gene_id"])
    
    df = df.sort_values(["chrom", "start"])
    
    # Save the selected rows to a new file
    
    return df.to_csv(sep="\t", header=False, index=False, doublequote=False, quoting=3, quotechar="",  escapechar="\\")


def analyze_target_ids(df, no_sgRNASet):
    # Count the occurrences of each transcript_id
    counts = df['target_id'].value_counts()

    # Create a DataFrame from the counts
    counts_df = counts.reset_index()
    counts_df.columns = ['target_id', 'count']


    counts = df['target_id'].value_counts()

    median_count = counts.median()
    min_count = counts.min()
    max_count = counts.max()

    print(f"\n\tMedian number of sgRNAs per target: {median_count}")
    print(f"\tMinimum number of sgRNAs per target: {min_count}")
    print(f"\tMaximum number of sgRNAs per target: {max_count}")
    print(f"\tNumber of targets with 0 sgRNA guides: {len(no_sgRNASet)}\n\n")

    return counts


    
def validate_scored_sgRNA_columns(args):
    required_columns = ['#chr', 'start', 'stop', 'sequence', 'pam', 'strand']

    filename = args.scoredGuides
    score_column=args.cleavage_column
    specificity_column = args.specificity_column
    
    with open(filename, 'r') as file:
        # Read the first line for the header
        header_line = file.readline().strip()
        # Split the line into columns based on tab delimiter
        columns = header_line.split('\t')

        # Check the fixed columns
        if score_column not in columns:
            raise ValueError("\t" + score_column + "not present")

        # Check the remaining columns are of the form "specificity..."
        if specificity_column is None:
            specificity_columns = [s for s in columns if "specificity" in s]
            if specificity_columns:
                specificity_column = specificity_columns[0]
            else:
                raise ValueError("\tNo off-target column specified and no column with 'specificity' found")
        elif specificity_column not in columns:
            raise ValueError("\t" + specificity_column + "not present")

    
    return columns, specificity_column, score_column 


def sgRNA_to_tscript(sgRNAs, mode, targets, header):
    
    df = pd.read_table(pyb.BedTool(sgRNAs).intersect(targets, wo=True).fn, sep="\t", header=None)
    df.columns = header + df.columns[len(header):].tolist()

    if mode == "transcript":
        df[['target_id', 'gene_id']] = df.iloc[:, len(df.columns) - 2].apply(extract_ids).apply(pd.Series)
    else:
        # assign an arbitrary number as it will be dropped later
        df[[999, 'target_id']] = df.iloc[:, len(df.columns) - 2].apply(extract_ids).apply(pd.Series)
    # Dropping integer-labeled columns
    int_columns = [col for col in df.columns if isinstance(col, int)]
    df = df.drop(int_columns, axis=1)
    #print(df.head())    

    df_without = pd.read_table(pyb.BedTool(targets).intersect(sgRNAs, v=True).fn, sep="\t", header=None)
    df_without[['transcript', 'gene']] = df_without.iloc[:, len(df_without.columns) - 1].apply(extract_ids).apply(pd.Series)
    
    set_without = set(df_without[mode])

    return df,set_without

def sgRNA_to_bed(sgRNAs, targets, header, target_header):
    
    df = pd.read_table(pyb.BedTool(sgRNAs).intersect(targets, wo=True).fn, sep="\t", header=None)

    #drop the last column (number of overlapping bases)
    df = df.drop(df.columns[-1], axis=1)

    num_columns = len(df.columns)

    sgRNAheader = header + target_header

    # Determine the number of columns to rename (minimum of the two lengths)
    num_rename = min(num_columns, len(sgRNAheader))

    # Create a dictionary mapping current column names to new names
    rename_dict = dict(zip(df.columns[:num_rename], sgRNAheader[:num_rename]))

    # Rename the columns
    df.rename(columns=rename_dict, inplace=True)

    df_without = pd.read_table(pyb.BedTool(targets).intersect(sgRNAs, v=True).fn, sep="\t", header=None)
    set_without = set(df_without[3])

    return df, set_without


def filter_df_by_column(df, column, min_value):
    print(f"\tFiltering by {column}, cut-off:   {min_value}")
    df = df[df[column] >= min_value]
    return df

def group_and_minimize(df, rank_column, num_to_keep):
    df = df.sort_values(by=['target_id', rank_column], ascending=[True, False])
    if num_to_keep == -1:
        print("\t--num_to_keep=-1, all sgRNAs returned")
        return df

    print(f"\tKeeping top {num_to_keep} sgRNAs per target, according to {rank_column}")
    # Group by target_id and keep the top entries as specified by number_of_guides
    df = df.groupby('target_id').head(num_to_keep)
    return df


def add_weighted_average_column(df, args):
    df['combined_weighted'] = ((df[args.specificity_column] * args.specificity_weight + df[args.cleavage_column] * args.cleavage_weight) / \
        (args.specificity_weight + args.cleavage_weight)).round(4)
    return df

def select_guides(df, rank_column, min_spacing):
    df =  df.sort_values(by=['target_id', rank_column], ascending=[True, False])

    if min_spacing == 0:
        print("\t--min_spacing=0, no filtering applied")
        return df

    print(f"\tEnforcing minimum spacing of {min_spacing} nucleotides between guides")

    def is_overlapping_or_close(new_start, new_stop, selected_guides):
        for start, stop in selected_guides:
            if (new_start <= stop + min_spacing) and (new_stop >= start - min_spacing):
                return True
        return False

    def select_from_group(group):
        selected_guides = []  # List to store tuples of (start, stop)

        for _, row in group.iterrows():
            if not is_overlapping_or_close(row['start'], row['stop'], selected_guides):
                selected_guides.append((row['start'], row['stop']))

        return group[group.apply(lambda row: (row['start'], row['stop']) in selected_guides, axis=1)]

    # Assuming df is already grouped by 'target_id'
    return df.groupby('target_id', group_keys=False).apply(select_from_group)

# might be a faster implementation for spacing guides.
# from intervaltree import Interval, IntervalTree
# def select_guides(df, min_spacing):
#     def select_from_group(group):
#         tree = IntervalTree()
#         selected_indices = []

#         for index, row in group.iterrows():
#             start, stop = row['start'], row['stop']
#             if not tree.overlaps(start - min_spacing, stop + min_spacing):
#                 tree.add(Interval(start, stop))
#                 selected_indices.append(index)

#         return group.loc[selected_indices]

#     return df.groupby('target_id', group_keys=False).apply(select_from_group)

def update_no_targets(original_target_ids, target_df):
    updated_target_ids = set(target_df['target_id'].unique())
    removed_target_ids = original_target_ids - updated_target_ids
    return removed_target_ids

def create_histogram(df, output_pdf_path):
    # Count the occurrences of each 'target_id'
    target_id_counts = df['target_id'].value_counts()

    # Count the frequencies of these occurrences
    frequency_counts = target_id_counts.value_counts().sort_index()

    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    frequency_counts.plot(kind='bar')
    plt.title('Frequency of Target ID Counts')
    plt.xlabel('Number of sgRNAs per Target ID')
    plt.ylabel('Unique target IDs with given number of sgRNAs')
    plt.tight_layout()  # Adjust layout to prevent clipping of labels

    # Save the figure to a PDF file
    plt.savefig(output_pdf_path)
    plt.close()


def save_notarget_set(noSgRNATargetSet, filename):
    
    sorted_strings = sorted(noSgRNATargetSet)
    with open(filename, "w") as file:
        for item in sorted_strings:
            file.write(item + "\n")

def make_output_dict(outDir, tempOutputDir, prefix):

    outputDict = {
        "df1" : tempOutputDir + prefix + "1_sgRNAsTargets_preFilter.tsv",
        "df2" : tempOutputDir + prefix + "2_sgRNAsTargets_specFilter.tsv",
        "df3" : tempOutputDir + prefix + "3_sgRNAsTargets_cleaveFilter.tsv",
        "df4" : tempOutputDir + prefix + "4_sgRNAsTargets_spacingFilter.tsv",
        "df5" : outDir + prefix + "FinalRankedsgRNAsTargets.tsv",    
        "empty" : outDir + prefix + "TargetsWith0sgRNAs.tsv",
        "counts" : outDir + prefix + "TargetsgRNACounts.tsv"
    }

    return outputDict

if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    args = parse_arguments()

    # make "tmp" directory
    outputDir = "./sgRNAs/RankedsgRNAs/"
    Path(outputDir).mkdir(parents=True, exist_ok=True)
    tempOutputDir = "./sgRNAs/RankedsgRNAs/tmp/"
    Path(tempOutputDir).mkdir(parents=True, exist_ok=True)


    # Load sgRNAs.bed file
    header, args.specificity_column, score_column = validate_scored_sgRNA_columns(args = args)
    sgRNAs = pyb.BedTool(args.scoredGuides)
    
    # Determine target file type
    targetFileType = args.targets.split(".")[-1].lower()
    targets = None
    targetsWithoutsgRNAs = None

    outputDict = make_output_dict(outDir=outputDir, tempOutputDir=tempOutputDir, prefix=args.prefix)
    if args.plot_histogram or args.output_all:
        args.output_all = True

    if targetFileType in ["gff", "gtf", "gff3", "gff2"]:
        print(f'\n\t{args.targets} is {targetFileType.upper()} format')
        trimmed = tempOutputDir + "trimmedAnnotation_" + args.targets
        with open(trimmed, 'w') as f:
            f.writelines(process_gtf(input_file = args.targets, feature = args.feature, percentiles = [0,args.percentile_range * 100]))
        targets = pyb.BedTool(trimmed)

        finalsgRNAs, targetsWithoutsgRNAs = sgRNA_to_tscript(sgRNAs=sgRNAs, mode=args.target_mode, targets=targets, header=header)

    else:
        if targetFileType not in ["bed"]:
            print(f'{args.targets} file extension unrecognized, for expected behavior pass GFF/GTF/BED file.\n \
            Defaulting to BED behavior')
        targets = pyb.BedTool(args.targets)
        target_header = ["target_chrom", "target_start", "target_stop", "target_id", 
                        "target_score", "target_strand", "target_big_start", "target_big_stop"]
        
        finalsgRNAs, targetsWithoutsgRNAs = sgRNA_to_bed(sgRNAs=sgRNAs, targets=targets, header=header, target_header=target_header)

    
    original_target_ids = set(finalsgRNAs['target_id'].unique())
    finalsgRNAs = add_weighted_average_column(df=finalsgRNAs, args=args)
    if args.output_all: finalsgRNAs.to_csv(outputDict["df1"], sep='\t', index=False)
    if args.plot_histogram: create_histogram(finalsgRNAs, outputDict["df1"][:-3] + "pdf")

    print("\n\n\tPrior to filtering: ")
    _ = analyze_target_ids(finalsgRNAs, targetsWithoutsgRNAs) 
    
    # Filter by specificity column
    finalsgRNAs = filter_df_by_column(df=finalsgRNAs, column=args.specificity_column, min_value=args.min_specificity)
    _ = analyze_target_ids(finalsgRNAs, targetsWithoutsgRNAs.union(update_no_targets(original_target_ids,finalsgRNAs)))
    if args.output_all: finalsgRNAs.to_csv(outputDict["df2"], sep='\t', index=False)
    if args.plot_histogram: create_histogram(finalsgRNAs, outputDict["df2"][:-3] + "pdf")
    
    # Filter by cleavage column
    finalsgRNAs = filter_df_by_column(df=finalsgRNAs, column=args.cleavage_column, min_value=args.min_cleavage)
    _ = analyze_target_ids(finalsgRNAs, targetsWithoutsgRNAs.union(update_no_targets(original_target_ids,finalsgRNAs)))
    if args.output_all: finalsgRNAs.to_csv(outputDict["df3"], sep='\t', index=False)
    if args.plot_histogram: create_histogram(finalsgRNAs, outputDict["df3"][:-3] + "pdf")


    # Filter by spacing
    finalsgRNAs = select_guides(df=finalsgRNAs, rank_column=args.rank_column, min_spacing=args.min_spacing)
    _ = analyze_target_ids(finalsgRNAs, targetsWithoutsgRNAs.union(update_no_targets(original_target_ids,finalsgRNAs)))
    if args.output_all: finalsgRNAs.to_csv(outputDict["df4"], sep='\t', index=False)
    if args.plot_histogram: create_histogram(finalsgRNAs, outputDict["df4"][:-3] + "pdf")


    # Filter by number_of_guides per target
    finalsgRNAs = group_and_minimize(finalsgRNAs, args.rank_column, args.number_of_guides)
    counts = analyze_target_ids(finalsgRNAs, targetsWithoutsgRNAs.union(update_no_targets(original_target_ids,finalsgRNAs)))
    finalsgRNAs.to_csv(outputDict["df5"], sep='\t', index=False)
    if args.plot_histogram: create_histogram(finalsgRNAs, outputDict["df5"][:-3] + "pdf")
    save_notarget_set(noSgRNATargetSet = targetsWithoutsgRNAs.union(update_no_targets(original_target_ids,finalsgRNAs)), filename=outputDict["empty"])
    counts.to_csv(outputDict["counts"], sep='\t', index=True)





    # TODO
    # add optional output at each step
    # optional plotting