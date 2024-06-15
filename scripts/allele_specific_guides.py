#!/usr/bin/env python3

import pandas as pd
import argparse

'''
    ./allele_specific_guides.py path/to/first_sgrna.bed path/to/second_sgrna.bed -d 1
'''

def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def find_unique_sgRNAs(sgrna_bed_1, sgrna_bed_2, min_distance=1):
    pd.set_option('display.max_columns', None)

    file1 = sgrna_bed_1.split("/")[-1].split(".")[0]
    file2 = sgrna_bed_2.split("/")[-1].split(".")[0]
    
    df1 = pd.read_csv(sgrna_bed_1, delimiter='\t')
    df2 = pd.read_csv(sgrna_bed_2, delimiter='\t')

    df1['id'] = df1['id,sequence,pam,chromosome,position,sense'].apply(lambda x: x.split(',')[0])
    df2['id'] = df2['id,sequence,pam,chromosome,position,sense'].apply(lambda x: x.split(',')[0])

    df1['sequence'] = df1['id,sequence,pam,chromosome,position,sense'].apply(lambda x: x.split(',')[1])
    df2['sequence'] = df2['id,sequence,pam,chromosome,position,sense'].apply(lambda x: x.split(',')[1])

    merged_df = df1.merge(df2, on="id", suffixes=('_' + file1, '_' + file2))

    unique_df1 = df1[~df1['id'].isin(merged_df['id'])].copy()
    unique_df2 = df2[~df2['id'].isin(merged_df['id'])].copy()

    unique_df1.drop(['id', 'sequence'], axis=1, inplace=True)
    unique_df2.drop(['id', 'sequence'], axis=1, inplace=True)

    # Output these unique entries to new files
    unique_df1.to_csv(file1 + '_specificPAMs.bed', sep='\t', index=False)
    unique_df2.to_csv(file2 + '_specificPAMs.bed', sep='\t', index=False)

    # Calculate Hamming distances
    merged_df['hamming_dist'] = merged_df.apply(
        lambda row: hamming_distance(row['sequence_' + file1], row['sequence_' + file2]), axis=1)
    merged_df['context_hamming_dist'] = merged_df.apply(
        lambda row: hamming_distance(row['context_' + file1], row['context_' + file2]), axis=1)

    filtered_df = merged_df[merged_df['hamming_dist'] >= min_distance]

    df1 = filtered_df[['#chr_'+file1, 'start_'+file1, 'stop_'+file1,'id,sequence,pam,chromosome,position,sense_'+file1,
                       'context_'+file1,'strand_'+file1, 'sequence_'+file2,'context_'+file2,'hamming_dist','context_hamming_dist']]
    df1_renamedict = {'id,sequence,pam,chromosome,position,sense_'+file1:'id,sequence,pam,chromosome,position,sense',
                      'context_'+file1:'context','#chr_'+file1:'#chr','start_'+file1:'start','stop_'+file1:'stop','strand_'+file1:'strand'}
    df1 = df1.rename(columns=df1_renamedict)

    df2 = filtered_df[['#chr_'+file2,'start_'+file2,'stop_'+file2,'id,sequence,pam,chromosome,position,sense_'+file2,
                    'context_'+file2,'strand_'+file2, 'sequence_'+file1,'context_'+file1,'hamming_dist','context_hamming_dist']]
    df2_renamedict = {'id,sequence,pam,chromosome,position,sense_'+file2:'id,sequence,pam,chromosome,position,sense',
                      'context_'+file2:'context','#chr_'+file2:'#chr','start_'+file2:'start','stop_'+file2:'stop','strand_'+file2:'strand'}
    df2 = df2.rename(columns=df2_renamedict)

    df1.to_csv( file1 + '_with_dist.bed', sep='\t', index=False)
    df2.to_csv( file2 + '_with_dist.bed', sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find unique sgRNAs and calculate hamming distances.')

    parser.add_argument('sgrna_bed_1', type=str, help='Path to the first sgRNA bed file.')
    parser.add_argument('sgrna_bed_2', type=str, help='Path to the second sgRNA bed file.')
    parser.add_argument('-d', '--min_distance', type=int, default=1, help='Minimum Hamming distance to filter sgRNAs. Default is 1.')

    args = parser.parse_args()

    find_unique_sgRNAs(args.sgrna_bed_1, args.sgrna_bed_2, args.min_distance)