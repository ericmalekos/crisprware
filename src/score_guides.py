#!/usr/bin/env python3

'''
    This script takes the output from FiltersgRNABED and scores it
    with Ruleset3 and Guidescan2. 

    ./crispomics/utils/ScoreGuides -k sgRNAs/NCAM_sgRNAs.bed -i chr21Index/ --tracr Chen2013 --threads 8
'''

import pandas as pd
import subprocess
import argparse
from rs3.seq import predict_seq
from scipy.stats import norm
from pathlib import Path
# import os
# import time

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
        default='Chen2013',
        required=True, 
        choices=['Hsu2013', 'Chen2013'], 
        help="TracrRNA version for cleavage scoring. \
            Either 'Hsu2013' or 'Chen2013', see https://github.com/gpp-rnd/rs3 \
            for details."
    )

    parser.add_argument(
        "--threads", 
        type=int, 
        default=8, 
        help="Number of threads [default: 8]"
    )
    
    parser.add_argument(
        "--mismatches", 
        type=int, 
        default=3, 
        help="Number of mismatches for Guidescan2 off-target scoring [default: 3]"
    )

    parser.add_argument(
        "-d", "--drop_duplicates", 
        help="Drop exact duplicate sgRNAs before scoring to save time.\
             Set flag to retain duplicates. [default: True]",
        default=True,
        action="store_false"
    )

    parser.add_argument(
        "--min_rs3", 
        type=float, 
        default=-2, 
        help="Minimum cleavage RS3 score. RS3 cleavage scores are formatted \
            as z-scores, so this is similar to a standard deviation cutoff. \
            [default: -2]"
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
    return parser.parse_args()

def guideScanScoring(guideCSV, output, guideIndex, threads = 2, mismatches = 3):
    
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
        '--mode',
        'succinct',
        '--kmers-file',
        guideCSV,
        '--output',
        output,
        guideIndex
    ]

    subprocess.run(cmd, check=True)

    # read the csv file
    print("\n\tReading in " + output)
    gscanDF = pd.read_csv(output)

    # drop duplicate rows and keep the first occurrence
    gscanDF = gscanDF.drop_duplicates(subset='id', keep='first')

    gscanDF['specificity'] = gscanDF['specificity'].round(3)

    db_name = Path(guideIndex).parts[-1]

    gscanDF.rename(columns={'specificity': 'specificity_' + db_name}, inplace=True)

    gscanDF = gscanDF.drop(['sequence', 'match_chrm',  'match_position', 'match_strand',  'match_distance'], axis = 1)

    return gscanDF

def cleavageScoring(sgRNADF, tracr = 'Chen2013', threads = 2, chunk_size = 2500000, minStdDev = None):

    print(f"\n\tBeginning RS3 cleavage scoring\n\tIf memory constrained reduce '--chunk_size'\n")

    # guideScores = predict_seq(sgRNAlist, sequence_tracr='Hsu2013', n_jobs=12)

    sgRNAlist = sgRNADF['context'].tolist() # 1709 seconds chr21, 12 cores 

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
    sgRNADF.loc[:, 'id'] = sgRNADF['id,sequence,pam,chromosome,position,sense'].str.split(',').str[0]


    return sgRNADF


def check_files_exist(index):
    """Check the existence of the three required files for a given index."""
    files = [f"{index}.reverse", f"{index}.forward", f"{index}.gs"]
    for file in files:
        if not Path(file).exists():
            raise FileNotFoundError(f"\n\n\tRequired file {file} not found for index {index}. \
                                    \n\tMake sure \n\t\t{index}.reverse \n\t\t{index}.forward \n\t\t{index}.gs \n\texist\n")


def main():
    
    args = parse_arguments()

    for index in args.guidescan2_indices:
        check_files_exist(index)

    # p = Path(sgRNA_output_path)
    # Path("./" + "/".join(p.parts[:-1])).mkdir(parents=True, exist_ok=True)

    output_path = "./sgRNAs/scoredSgRNAs/" + args.output_prefix
    tmp_path = "./sgRNAs/scoredSgRNAs/tmp/"  + args.output_prefix
    o, t = Path(output_path), Path(tmp_path)
    #Path("./" + "/".join(o.parts[:-1])).mkdir(parents=True, exist_ok=True)
    Path("./" + "/".join(t.parts[:-1])).mkdir(parents=True, exist_ok=True)

    sgRNADF = pd.read_csv(args.sgrna_bed, delimiter='\t', header=0 )

    if args.drop_duplicates:
        print(f'\n\tBefore dropping duplicates:\t{len(sgRNADF)}')
        sgRNADF['sgRNA'] = sgRNADF.iloc[:, 3].str.split(',').str[1]

        # Step 3 & 4: Find unique and duplicate k-mers, then filter the DataFrame to keep only the unique ones
        unique_mask = ~sgRNADF['sgRNA'].duplicated(keep=False)
        sgRNADF = sgRNADF[unique_mask]

        print(f'\tAfter dropping duplicates:\t{len(sgRNADF)}\n')

    sgRNADF = cleavageScoring(sgRNADF = sgRNADF, tracr = args.tracr, chunk_size=args.chunk_size, threads = args.threads, minStdDev = args.min_rs3)
    
    print(f'\n\tAfter dropping RS3 cleavage scores below {args.min_rs3}:\t{len(sgRNADF)}\n')

    num_chunks = len(sgRNADF) // args.chunk_size + (1 if len(sgRNADF) % args.chunk_size else 0)
    guidescan_dfs = []

    for gscanIndex in args.guidescan2_indices:
        print(f"\n\tBeginning Guidescan2 specificity scoring against " + gscanIndex + "\n\tIf memory constrained reduce '--chunk_size'\n")
        guidescan_chunk_dfs = []
        # chunk up the input and save it in a form compatible with guidescan processing
        for i, (_, chunk) in enumerate(sgRNADF.groupby(sgRNADF.index // args.chunk_size)):

            input = tmp_path + f'gscanInput.{i + 1}.csv'
            output = input.replace("Input", "Output")
            print(input)
            chunk[['id,sequence,pam,chromosome,position,sense']].to_csv("./" + input, sep='\t', index=False)
            #guidescan_files.append(tuple(input, output))
            guidescan_chunk_dfs.append(guideScanScoring(guideCSV = input, output = output, 
                                                            guideIndex = gscanIndex, 
                                                            threads = args.threads, 
                                                            mismatches = args.mismatches))
        
        # rebuild the initial df by concatenating the chunks
        guidescan_dfs.append(pd.concat(guidescan_chunk_dfs, ignore_index=True))


    for df in guidescan_dfs:
        sgRNADF = sgRNADF.merge(df, on='id', how='inner')
    
    sgRNADF[['sequence', 'pam']] = sgRNADF['id,sequence,pam,chromosome,position,sense'].str.split(',', expand=True).iloc[:, 1:3]

    sgRNADF = sgRNADF.drop(['id,sequence,pam,chromosome,position,sense', 'id'], axis = 1)
    
    desired_cols = ['#chr', 'start', 'stop', 'sequence', 'pam', 'strand', 'context',
                   'rs3_z_score', 'rs3_cdf','rs3_score_norm']
    specificity_cols = [col for col in sgRNADF.columns if "specificity" in col]

    sgRNADF = sgRNADF[desired_cols + specificity_cols]

    sgRNADF.to_csv(output_path + "AllScoredGuides.tsv", sep="\t", index=False)


    ### TODO
    # Add a flag that declines rs3 scoring
    #   in that case need to rethink the "desired columns" variable

if __name__ == "__main__":
    main()