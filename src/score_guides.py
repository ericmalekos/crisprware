#!/usr/bin/env python3

'''
    This script takes the output from generate_guides and adds
    with  on-target Ruleset3 and off-target Guidescan2 scores.

    score_guides --grna_bed tests/test_data/gRNA_test.bed \
        --guidescan2_indices tests/test_data/Hg38_chr21_gscan_index/chr21.fa.index\
        --tracr Chen2013 --threads 2 --output_prefix tests/test_output/Hg38_chr21

'''

import pandas as pd
import subprocess
import argparse
from rs3.seq import predict_seq
#from scipy.stats import norm
from pathlib import Path
from os import remove
from utils.utility_functions import create_output
from utils.dna_sequence_functions import map_ambiguous_sequence

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return x

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Scores gRNAs from generate_guides."
    )

    parser.add_argument(
        "-b", "--grna_bed",
        type=str,
        help="grnas.bed ouput of generate_guides.",
        required=True
    )

    parser.add_argument(
        "-i", "--guidescan2_indices",
        type=str,
        help="One or more, space-separate Guidescan2 indices. \
            A specificity score will be calculated against each index separately.",
        required=False,
        nargs='*'
    )

    parser.add_argument(
        "--tracr",
        type=str,
        default=None,
        required=False,
        choices=['Hsu2013', 'Chen2013', 'both'],
        help="TracrRNA version for cleavage scoring. \
            Either 'Hsu2013' or 'Chen2013' or 'both', see https://github.com/gpp-rnd/rs3 \
            for details."
    )

    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of threads [default: 8]"
    )

    parser.add_argument(
        "--threshold",
        type=int,
        default=2,
        help="Threshold for Guidescan2 off-target hits. If off-targets are found this distance away \
         the gRNA will be discarded, i.e. set to 2 to discard any guides with a 0, 1 or 2 mismatches \
         from another PAM adjacent sequence. --threshold=-1 to retain all guides [default: 2]"
    )

    parser.add_argument(
        "--mismatches",
        type=int,
        default=3,
        help="Number of mismatches for Guidescan2 off-target scoring [default: 3]"
    )

    parser.add_argument(
        "--rna_bulges",
        type=int,
        default=0,
        help="RNA bulges for Guidescan2 off-target scoring [default: 0]"
    )

    parser.add_argument(
        "--dna_bulges",
        type=int,
        default=0,
        help="DNA bulges for Guidescan2 off-target scoring [default: 0]"
    )

    parser.add_argument(
        "--mode",
        type=str,
        default='succinct',
        choices=['succinct', 'complete'],
        help="Whether Guidescan2 temporary output should be succinct or complete mode [default: 0]"
    )

    parser.add_argument(
        "--alt_pams",
        type=str,
        help="One or more, space-separate alternative pams for off-target \
            consideration. e.g. NAG",
        required=False,
        nargs='*'
    )

    parser.add_argument(
        "-d", "--drop_duplicates",
        help="Drop exact duplicate gRNAs before scoring to save time.\
             Set flag to retain duplicates. [default: True]",
        default=True,
        action="store_false"
    )

    parser.add_argument(
        "--skip_rs3",
        help="Set flag to skip RS3 scoring [default: False]",
        default=False,
        action="store_true"
    )

    parser.add_argument(
        "--skip_gs2",
        help="Set flag to skip Guidescan2 scoring [default: False]",
        default=False,
        action="store_true"
    )

    parser.add_argument(
        "--min_rs3",
        type=float,
        default=None,
        help="Minimum cleavage RS3 score. RS3 cleavage scores are formatted \
            as z-scores, so this is interpreted as a standard deviation cutoff. \
            Functionality also available in rank_guides.py. Applying at this \
            stage can increase speed by filtering before off-target scoring. \
            [default: None]"
    )

    parser.add_argument(
        "--chunk_size",
        type=int,
        default=100000,
        help="Number of gRNAs to hold in memory for cleavage scoring and off-target filtering. \
        Reduce if memory constrained. Increasing may improve runtime [default: 100000]"
    )

    parser.add_argument(
        "-k", "--keep_tmp",
        help="Set flag to keep temporary Guidescan2 output [default: False]",
        default=False,
        action="store_true"
    )

    parser.add_argument(
        "-o", "--output_directory",
        help="Path to output. [default: current directory]",
        type=str,
        default=""
    )

    return parser.parse_args()




def gscan_scoring(guideCSV, output, guideIndex, \
                    threads = 2, mismatches = 3, rna_bulges = 0, dna_bulges = 0, \
                    threshold = 2, mode = "succinct", alt_pam = "", keep_gscan = False):
    """
    Executes GuideScan's enumerate command with specified parameters and processes the output.

    Parameters:
    - guideCSV (str): The path to the input CSV file containing guide RNA sequences.
    - output (str): The path where the output CSV file will be saved.
    - guideIndex (str): The index used for off-target scoring.
    - threads (int, optional): The number of threads to use for the GuideScan command. Defaults to 2.
    - mismatches (int, optional): The maximum number of mismatches allowed. Defaults to 3.
    - rna_bulges (int, optional): The number of RNA bulges allowed. Defaults to 0.
    - dna_bulges (int, optional): The number of DNA bulges allowed. Defaults to 0.
    - threshold (int, optional): The threshold parameter for removing low mismatches. Defaults to 2.
    - mode (str, optional): The mode of operation for the GuideScan command. Defaults to "succinct".
    - alt_pam (str, optional): An alternative PAM sequence to be considered. Defaults to an empty string.

    Returns:
    - pandas.DataFrame: A DataFrame containing the processed output from the GuideScan command.
    """

    cmd = [
        'guidescan',
        'enumerate',
        '--max-off-targets',
        '-1',
        '--threads',
        str(threads),
        '--mismatches',
        str(mismatches),
        '--format',
        'csv',
        '--rna-bulges',
        str(rna_bulges),
        '--dna-bulges',
        str(dna_bulges),
        '--threshold',
        str(threshold),
        '--mode',
        mode,
        '--alt-pam',
        'None',
        '--kmers-file',
        guideCSV,
        '--output',
        output,
        guideIndex
    ]

    subprocess.run(cmd, check=True)

    # read the csv file
    gscanDF = pd.read_csv(output)

    #gscanDF = process_gscan(gscanDF, guideIndex)
    if gscanDF.empty:
        print(f'\n\tGuidescan2 output is empty, consider lowering "--threshold"\n')

    gscanDF = gscanDF.drop_duplicates(subset='id', keep='first')

    gscanDF['specificity'] = gscanDF['specificity'].round(4)

    db_name = Path(guideIndex).parts[-1]

    gscanDF.rename(columns={'specificity': 'specificity_' + db_name}, inplace=True)

    gscanDF = gscanDF.drop(['sequence', 'match_chrm',  'match_position', 'match_strand', 'match_distance'], axis = 1)

    if not keep_gscan:
        remove(output)
        remove(guideCSV)

    return gscanDF

def compute_rs3_scores(gRNAlist, tracr, threads, chunk_size):
    """
    Compute RS3 scores for a list of gRNAs with the specified tracrRNA sequence.

    Parameters:
    - gRNAlist (list): List of gRNA 30-nt contexts.
    - tracr (str): The tracrRNA sequence identifier to be used in cleavage scoring.
    - threads (int): The number of parallel jobs to run for cleavage scoring.
    - chunk_size (int): The size of chunks to split the gRNA list into for processing.

    Returns:
    - list: List of RS3 scores for the input gRNA sequences.
    """
    gRNAScores = []

    # Iterate over big_list in chunks of size chunk_size
    for i in range(0, len(gRNAlist), chunk_size):
        sublist = gRNAlist[i:i + chunk_size]
        processed_sublist = predict_seq(sublist, sequence_tracr=tracr, n_jobs=threads)
        gRNAScores.extend(processed_sublist)

    return gRNAScores

def cleavage_scoring(gRNADF, tracr, threads = 2, chunk_size = 200000, minStdDev = None):
    """
    Computes RS3 cleavage scores for gRNAs using the specified tracrRNA sequence and filters based on minimum standard deviation.

    Parameters:
    - gRNADF (pandas.DataFrame): DataFrame containing gRNA sequences with a column named 'context' for gRNA 30-nt contexts.
    - tracr (str): The tracrRNA sequence identifier to be used in cleavage scoring.
    - threads (int, optional): The number of parallel jobs to run for cleavage scoring. Defaults to 2.
    - chunk_size (int, optional): The size of chunks to split the gRNA list into for processing, to manage memory usage. Defaults to 2000000.
    - minStdDev (float or None, optional): The minimum standard deviation threshold for filtering gRNAs based on their RS3 z-score.
      gRNAs with a z-score below this threshold will be excluded. Defaults to None, which disables filtering.

    Returns:
    - pandas.DataFrame: The input DataFrame augmented with 'RS3_score' and 'rs3_cdf' columns, and optionally filtered
      based on the 'RS3_score' threshold.
    """

    print(f"\n\tBeginning RS3 cleavage scoring\n\tIf memory constrained reduce '--chunk_size'\n")

    gRNAlist = gRNADF['context'].tolist()

    if tracr == 'both':
        gRNAScores_Hsu2013 = compute_rs3_scores(gRNAlist, "Hsu2013", threads, chunk_size)
        gRNAScores_Chen2013 = compute_rs3_scores(gRNAlist, "Chen2013", threads, chunk_size)

        gRNADF['RS3_score_Hsu2013'] = gRNAScores_Hsu2013
        gRNADF['RS3_score_Hsu2013'] = gRNADF['RS3_score_Hsu2013'].round(4)
        # gRNADF['rs3_cdf_Hsu2013'] = norm.cdf(gRNADF['RS3_score_Hsu2013'])
        # gRNADF['rs3_cdf_Hsu2013'] = gRNADF['rs3_cdf_Hsu2013'].round(4)

        gRNADF['RS3_score_Chen2013'] = gRNAScores_Chen2013
        gRNADF['RS3_score_Chen2013'] = gRNADF['RS3_score_Chen2013'].round(4)
        # gRNADF['rs3_cdf_Chen2013'] = norm.cdf(gRNADF['RS3_score_Chen2013'])
        # gRNADF['rs3_cdf_Chen2013'] = gRNADF['rs3_cdf_Chen2013'].round(4)
    elif tracr in ["Hsu2013", "Chen2013"]:
        gRNAScores = compute_rs3_scores(gRNAlist, tracr, threads, chunk_size)

        gRNADF['RS3_score_' + tracr] = gRNAScores
        gRNADF['RS3_score_' + tracr] = gRNADF['RS3_score_' + tracr].round(4)
        # gRNADF['rs3_cdf'] = norm.cdf(gRNADF['RS3_score'])
        # gRNADF['rs3_cdf'] = gRNADF['rs3_cdf'].round(4)

        if minStdDev:
            gRNADF = gRNADF[gRNADF['RS3_score_' + tracr] > minStdDev]

    gRNADF = gRNADF.copy()

    return gRNADF

def check_files_exist(index):
    """Check the existence of the three required files for a given index."""

    files = [f"{index}.reverse", f"{index}.forward", f"{index}.gs"]
    for file in files:
        if not Path(file).exists():
            raise FileNotFoundError(f"\n\n\tRequired file {file} not found for index {index}. \
                                    \n\tMake sure \n\t\t{index}.reverse \n\t\t{index}.forward \n\t\t{index}.gs \n\texist\n")

def get_alt_pams(pams):
    """
    Generates a string of alternative PAM sequences by expanding each ambiguous PAM sequence in the input list.

    Parameters:
    - pams (list of str): A list of PAM sequences. These sequences can include ambiguous nucleotide symbols
      that represent multiple possible nucleotides at a particular position.

    Returns:
    - str: A string containing all unique alternative PAM sequences derived from the input list, space-separated.
    """

    pamlist = []
    for apam in pams:
        pamlist += map_ambiguous_sequence(apam)
    pamstr = " ".join(set(pamlist))
    return pamstr

def main():

    args = parse_arguments()

    #print("ALT PAMS" + args.alt_pams)
    if not args.skip_gs2:
        for index in args.guidescan2_indices:
            check_files_exist(index)

    if not args.skip_rs3 and not args.tracr:
        raise ValueError("\n\tThe '--tracr' argument is required for RS3 scoring but was not provided.\n")

    pams = ""
    if args.alt_pams:
        pams=get_alt_pams(args.alt_pams)
        print(f'\n\tConsidering {pams} for off-target scoring \n')
    else:
        args.alt_pams = None

    #gRNA_output_path = "./" + args.output_prefix + "ScoredSgRNAs/" + args.output_prefix.split("/")[-1] + "ScoredSgRNAs.tsv"
    #tmp_path = create_output_directory(base_dir="./" + args.output_prefix + "ScoredSgRNAs/",output_prefix="tmp/")

    gRNA_output_path, tmp_path = create_output(args.grna_bed, outdir=args.output_directory, extension="scoredgRNA", stripped="_gRNA", tmp=True)
    gRNA_output_path += ".bed"
    #print(gRNA_output_path)
    #print(tmp_path)

    specificity_cols=[]

    gRNADF = pd.read_csv(args.grna_bed, delimiter='\t', header=0)
    final_columns = gRNADF.columns.tolist() + ['sequence']
    final_columns.remove('id,sequence,pam,chromosome,position,sense')

    if args.drop_duplicates:
        print(f'\n\tBefore dropping duplicates:\t{len(gRNADF)}')
        gRNADF['gRNA'] = gRNADF.iloc[:, 3].str.split(',').str[1]

        # Step 3 & 4: Find unique and duplicate k-mers, then filter the DataFrame to keep only the unique ones
        unique_mask = ~gRNADF['gRNA'].duplicated(keep=False)
        gRNADF = gRNADF[unique_mask]
        print(f'\tAfter dropping duplicates:\t{len(gRNADF)}\n')

    if not args.skip_rs3:
        if args.tracr == 'both':
            final_columns += ['RS3_score_Hsu2013', 'RS3_score_Chen2013']
        else:
            final_columns += ['RS3_score_' + args.tracr]
        gRNADF = cleavage_scoring(gRNADF = gRNADF,
                                   tracr = args.tracr,
                                   chunk_size=args.chunk_size,
                                   threads = args.threads,
                                   minStdDev = args.min_rs3)
        print(f'\n\tAfter dropping RS3 cleavage scores below {args.min_rs3}:\t{len(gRNADF)}\n')

    gRNADF.loc[:, 'id'] = gRNADF['id,sequence,pam,chromosome,position,sense'].str.split(',').str[0]
    guidescan_dfs = []

    if not args.skip_gs2:
        for gscanIndex in args.guidescan2_indices:
            print(f"\n\tBeginning Guidescan2 specificity scoring against " + gscanIndex + "\n\tIf memory constrained reduce '--chunk_size'\n")
            guidescan_chunk_dfs = []
            suffix_index = gscanIndex.split("/")[-1]
            # chunk up the input and save it in a form compatible with guidescan processing
            for i, (_, chunk) in enumerate(gRNADF.groupby(gRNADF.index // args.chunk_size)):

                input = tmp_path + f'{suffix_index}Input.{i + 1}.csv'
                print("input:"  + input)
                output = input.replace("Input", "Output")
                print(f"\n\tSaved Guidescan input file to {input}\n")
                chunk[['id,sequence,pam,chromosome,position,sense']].to_csv(input, sep='\t', index=False)
                guidescan_chunk_dfs.append(gscan_scoring(guideCSV=input, output=output,
                                                                guideIndex=gscanIndex,
                                                                threads=args.threads,
                                                                mismatches=args.mismatches,
                                                                rna_bulges=args.rna_bulges,
                                                                dna_bulges=args.dna_bulges,
                                                                mode=args.mode,
                                                                threshold=args.threshold,
                                                                alt_pam=pams,
                                                                keep_gscan=args.keep_tmp))

            guidescan_dfs.append(pd.concat(guidescan_chunk_dfs, ignore_index=True))

        for df in guidescan_dfs:
            gRNADF = gRNADF.merge(df, on='id', how='outer')
            #gRNADF = gRNADF.merge(df, on='id', how='inner')

        # If any rows have NaN for all specificity columns, drop those rows
        specificity_cols = [col for col in gRNADF.columns if "specificity" in col]
        rows_to_drop = gRNADF[specificity_cols].isna().all(axis=1)
        gRNADF = gRNADF.loc[~rows_to_drop]

    gRNADF[['sequence']] = gRNADF['id,sequence,pam,chromosome,position,sense'].str.split(',', expand=True).iloc[:, 1:2]

    gRNADF = gRNADF.drop(['id,sequence,pam,chromosome,position,sense', 'id'], axis = 1)

    gRNADF = gRNADF[final_columns + specificity_cols]

    gRNADF.to_csv(gRNA_output_path, na_rep='-1', sep="\t", index=False)
    print(f"\n\tSaved output file to {gRNA_output_path}\n")

if __name__ == "__main__":
    main()
