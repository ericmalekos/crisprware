#!/usr/bin/env python3

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write
import argparse

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def process_bed_file(input_file):
    # Read the BED file into a DataFrame
    df = pd.read_csv(input_file, sep='\t', comment='#', header=0, names=[
        'chr', 'start', 'stop', 'id', 'context', 'strand'])

    # Drop the fourth column
    df = df.drop(columns=['id'])

    # Reverse complement the context column if the strand is negative
    df.loc[df['strand'] == '-', 'context'] = df[df['strand'] == '-']['context'].apply(reverse_complement)
    
    df['start'] += 1
    #df.loc[df['strand'] == '-', 'start'] -= 1
    #df.loc[df['strand'] == '-', 'stop'] -= 1

    # Sort the DataFrame by chromosome and start position
    df = df.sort_values(by=['chr', 'start'])
    
    return df

def stitch_sequences(df):
    stitched_sequences = []
    current_chr = None
    current_start = current_end = None
    current_sequence = []

    for _, row in df.iterrows():
        chr, start, stop, context, strand = row

        if current_chr == chr and start <= current_end:
            # If the entries overlap, extend the end coordinate and concatenate only the novel portion of the sequence
            overlap_length = current_end - start + 1
            novel_context = context[overlap_length:]
            current_end = max(current_end, stop)
            current_sequence.append(novel_context)
        else:
            # If the entries do not overlap, save the current sequence and start a new one
            if current_sequence:
                stitched_sequences.append((current_chr, current_start, current_end, ''.join(current_sequence)))
            
            current_chr = chr
            current_start = start
            current_end = stop
            current_sequence = [context]
    
    # Append the last sequence
    if current_sequence:
        stitched_sequences.append((current_chr, current_start, current_end, ''.join(current_sequence)))
    
    return stitched_sequences

def create_fasta(stitched_sequences, output_file, reverse_flag):
    records = []
    for chr, start, end, sequence in stitched_sequences:
        if reverse_flag:
            sequence = sequence[::-1]
        header = f"{chr}_{start}_{end}"
        records.append(SeqRecord(Seq(sequence), id=header, description=""))
    
    with open(output_file, 'w') as fasta_out:
        write(records, fasta_out, 'fasta')
def main():
    parser = argparse.ArgumentParser(description='Process and stitch BED file sequences, then output to FASTA format.')
    parser.add_argument('input_file', type=str, help='Path to the input BED file')
    parser.add_argument('output_file', type=str, help='Path to the output FASTA file')
    parser.add_argument('--reverse', action='store_true', default=False, help='Reverse the sequence before outputting (default: False)')

    args = parser.parse_args()
    
    df = process_bed_file(args.input_file)
    stitched_sequences = stitch_sequences(df)
    create_fasta(stitched_sequences, args.output_file, args.reverse)
    print("Processing complete. The output has been saved to", args.output_file)

if __name__ == "__main__":
    main()
