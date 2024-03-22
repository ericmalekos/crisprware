#!/usr/bin/env python3

'''
    This script takes the output from generate_guides and adds
    with  on-target Ruleset3 and off-target Guidescan2 scores. 

    score_guides --sgrna_bed tests/test_data/sgRNA_test.bed \
        --guidescan2_indices tests/test_data/Hg38_chr21_gscan_index/chr21.fa.index\
        --tracr Chen2013 --threads 2 --output_prefix tests/test_output/Hg38_chr21

'''

import pandas as pd
import subprocess
import argparse
from rs3.seq import predict_seq
from scipy.stats import norm
from pathlib import Path
from os import remove
from utils.utility_functions import create_output_directory
from utils.dna_sequence_functions import map_ambiguous_sequence

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 1.0]")
    return x

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Scores sgRNAs from GenerateGuides."
    )

    parser.add_argument(
        "-b", "--sgrna_bed",
        type=str,
        help="sgrnas.bed ouput of GenerateGuides.",
        required=True
    )

    parser.add_argument(
        "-i", "--guidescan2_indices",
        type=str,
        help="One or more, space-separate Guidescan2 indices. \
            A specificity score will be calculated against each index separately.",
        required=True,
        nargs='*'
    )

    parser.add_argument(
        "--tracr", 
        type=str,
        default=None,
        required=False, 
        choices=['Hsu2013', 'Chen2013'], 
        help="TracrRNA version for cleavage scoring. \
            Either 'Hsu2013' or 'Chen2013', see https://github.com/gpp-rnd/rs3 \
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
         the sgRNA will be discarded, i.e. set to 2 to discard any guides with a 0, 1 or 2 mismatches \
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
        help="Drop exact duplicate sgRNAs before scoring to save time.\
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
        help="Number of sgRNAs to hold in memory for cleavage scoring and off-target filtering. \
        Reduce if memory constrained. Increasing may improve runtime [default: 100000]"
    )

    parser.add_argument(
        "-o", "--output_prefix",
        help="Prefix to use for files.",
        type=str,
        default="scoredOut_"
    )

    parser.add_argument(
        "-k", "--keep_tmp", 
        help="Set flag to keep temporary Guidescan2 output [default: False]",
        default=False,
        action="store_true"
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


def cleavage_scoring(sgRNADF, tracr, threads = 2, chunk_size = 200000, minStdDev = None):
    """
    Computes RS3 cleavage scores for sgRNAs using the specified tracrRNA sequence and filters based on minimum standard deviation.

    Parameters:
    - sgRNADF (pandas.DataFrame): DataFrame containing sgRNA sequences with a column named 'context' for sgRNA 30-nt contexts.
    - tracr (str): The tracrRNA sequence identifier to be used in cleavage scoring.
    - threads (int, optional): The number of parallel jobs to run for cleavage scoring. Defaults to 2.
    - chunk_size (int, optional): The size of chunks to split the sgRNA list into for processing, to manage memory usage. Defaults to 2000000.
    - minStdDev (float or None, optional): The minimum standard deviation threshold for filtering sgRNAs based on their RS3 z-score.
      sgRNAs with a z-score below this threshold will be excluded. Defaults to None, which disables filtering.

    Returns:
    - pandas.DataFrame: The input DataFrame augmented with 'rs3_z_score' and 'rs3_cdf' columns, and optionally filtered
      based on the 'rs3_z_score' threshold.
    """

    print(f"\n\tBeginning RS3 cleavage scoring\n\tIf memory constrained reduce '--chunk_size'\n")

    sgRNAlist = sgRNADF['context'].tolist()

    # process the list in chunks to reduce memory footprint
    sgRNAScores = []

    # Iterate over big_list in chunks of size chunk_size
    for i in range(0, len(sgRNAlist), chunk_size):
        sublist = sgRNAlist[i:i + chunk_size]
        processed_sublist = predict_seq(sublist, sequence_tracr=tracr, n_jobs=threads)
        sgRNAScores.extend(processed_sublist)

    sgRNADF['rs3_z_score'] = sgRNAScores
    sgRNADF['rs3_z_score'] = sgRNADF['rs3_z_score'].round(4)
    sgRNADF['rs3_cdf'] = norm.cdf(sgRNADF['rs3_z_score'])
    sgRNADF['rs3_cdf'] = sgRNADF['rs3_cdf'].round(4)

    if minStdDev:
        sgRNADF = sgRNADF[sgRNADF['rs3_z_score'] > minStdDev]

    sgRNADF = sgRNADF.copy()  

    return sgRNADF

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

    sgRNA_output_path = "./" + args.output_prefix + "ScoredSgRNAs/" + args.output_prefix.split("/")[-1] + "ScoredSgRNAs.tsv"
    tmp_path = create_output_directory(base_dir="./" + args.output_prefix + "ScoredSgRNAs/",output_prefix="tmp/")
    

    sgRNADF = pd.read_csv(args.sgrna_bed, delimiter='\t', header=0)
    final_columns = sgRNADF.columns.tolist() + ['sequence']
    final_columns.remove('id,sequence,pam,chromosome,position,sense')

    if args.drop_duplicates:
        print(f'\n\tBefore dropping duplicates:\t{len(sgRNADF)}')
        sgRNADF['sgRNA'] = sgRNADF.iloc[:, 3].str.split(',').str[1]

        # Step 3 & 4: Find unique and duplicate k-mers, then filter the DataFrame to keep only the unique ones
        unique_mask = ~sgRNADF['sgRNA'].duplicated(keep=False)
        sgRNADF = sgRNADF[unique_mask]
        print(f'\tAfter dropping duplicates:\t{len(sgRNADF)}\n')

    if not args.skip_rs3:
        final_columns += ['rs3_z_score', 'rs3_cdf']
        sgRNADF = cleavage_scoring(sgRNADF = sgRNADF,
                                   tracr = args.tracr,
                                   chunk_size=args.chunk_size,
                                   threads = args.threads,
                                   minStdDev = args.min_rs3)
        print(f'\n\tAfter dropping RS3 cleavage scores below {args.min_rs3}:\t{len(sgRNADF)}\n')

    sgRNADF.loc[:, 'id'] = sgRNADF['id,sequence,pam,chromosome,position,sense'].str.split(',').str[0]    
    guidescan_dfs = []

    for gscanIndex in args.guidescan2_indices:
        print(f"\n\tBeginning Guidescan2 specificity scoring against " + gscanIndex + "\n\tIf memory constrained reduce '--chunk_size'\n")
        guidescan_chunk_dfs = []
        # chunk up the input and save it in a form compatible with guidescan processing
        for i, (_, chunk) in enumerate(sgRNADF.groupby(sgRNADF.index // args.chunk_size)):

            input = tmp_path + f'gscanInput.{i + 1}.csv'
            output = input.replace("Input", "Output")
            print(f"\n\tSaved Guidescan input file to {input}\n")
            chunk[['id,sequence,pam,chromosome,position,sense']].to_csv("./" + input, sep='\t', index=False)
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
        sgRNADF = sgRNADF.merge(df, on='id', how='outer')
        #sgRNADF = sgRNADF.merge(df, on='id', how='inner')

    # If any rows have NaN for all specificity columns, drop those rows
    specificity_cols = [col for col in sgRNADF.columns if "specificity" in col]
    rows_to_drop = sgRNADF[specificity_cols].isna().all(axis=1)
    sgRNADF = sgRNADF.loc[~rows_to_drop]

    sgRNADF[['sequence']] = sgRNADF['id,sequence,pam,chromosome,position,sense'].str.split(',', expand=True).iloc[:, 1:2]

    sgRNADF = sgRNADF.drop(['id,sequence,pam,chromosome,position,sense', 'id'], axis = 1)

    sgRNADF = sgRNADF[final_columns + specificity_cols]

    sgRNADF.to_csv(sgRNA_output_path, na_rep='-1', sep="\t", index=False)
    print(f"\n\tSaved output file to {sgRNA_output_path}\n")

if __name__ == "__main__":
    main()