#!/usr/bin/env python3

'''
Usage:
./replace_snps.py ./chr19_GRCm39.primary_assembly.genome.fa output.fasta SNP_file1.txt SNP_file2.txt ...

'''

import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq

def preprocess_all_snps(snp_files):
    all_snps_dict = {}

    for snp_file_path in snp_files:
        snp_df = pd.read_csv(snp_file_path, sep='\t')

        for index, row in snp_df.iterrows():
            snp_chromosome = str(row['Chromosome'])

            if not snp_chromosome.startswith("chr"):
                snp_chromosome = "chr" + snp_chromosome

            if pd.isna(row['Coordinate (GRCm39)']):
                print(f"Warning: NaN value found for coordinate in {snp_file_path} at index {index}. Skipping this entry.")
                continue

            pos = int(row['Coordinate (GRCm39)']) - 1
            nucleotide = row['CAST/EiJ']

            if nucleotide not in {'A', 'T', 'C', 'G'}:
                print(f"Warning: Unexpected value '{nucleotide}' at position {pos} in {snp_file_path}. Skipping this entry.")
                continue

            key = (snp_chromosome, pos)
            all_snps_dict[key] = nucleotide

    return all_snps_dict

def replace_snps_in_fasta(fasta_path, output_path, snps_dict):
    # Initialize a flag to determine if we're processing the first chromosome
    first_chromosome = True

    # Read the FASTA file using Biopython
    for record in SeqIO.parse(fasta_path, "fasta"):
        mutable_seq = MutableSeq(str(record.seq))
        current_chromosome = record.id

        for (chromosome, pos), nucleotide in snps_dict.items():
            if chromosome == current_chromosome:
                try:
                    mutable_seq[pos] = nucleotide
                except Exception as e:
                    print(f"An error occurred while setting position {pos} to '{nucleotide}' in chromosome {chromosome}: {e}")

        record.seq = Seq(mutable_seq)

        # Determine the file mode based on whether it's the first chromosome or not
        mode = "w" if first_chromosome else "a"
        with open(output_path, mode) as output_handle:
            SeqIO.write(record, output_handle, "fasta")

        # Set the flag to False after processing the first chromosome
        first_chromosome = False

# Ensure there are at least 3 arguments (script name, fasta path, output path, and 1 or more SNP files)
if len(sys.argv) < 4:
    print("Usage: python replace_snps.py <path_to_fasta> <path_to_output_fasta> <path_to_snp_file1> [<path_to_snp_file2> ...]")
    sys.exit(1)

# Paths to the input and output files
fasta_path = sys.argv[1]
output_path = sys.argv[2]
snp_file_paths = sys.argv[3:]

# Run the function
all_snps_dict = preprocess_all_snps(snp_file_paths)
replace_snps_in_fasta(fasta_path, output_path, all_snps_dict)
