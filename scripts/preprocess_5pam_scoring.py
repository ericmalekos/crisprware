#!/usr/bin/env python3

import argparse

def reverse_sequence(seq):
    return seq[::-1]

def process_pam(pam):
    # Reverse the PAM and replace any non-ATCG letter with N
    return ''.join(['N' if base not in 'ATCG' else base for base in reverse_sequence(pam)])

def process_line(line, reverse_flag):
    columns = line.strip().split('\t')
    
    # Check if the line has the correct number of columns
    if len(columns) < 6:
        return None
    
    sequence_info = columns[3].split(',')
    sequence = sequence_info[1]
    #pam = sequence_info[2]
    pam=""
    sequence_info[3] = ""
    sequence_info[4] = ""
    sequence_info[5] = ""
    
    # Reverse the sequence if the reverse flag is set
    if reverse_flag:
        reversed_sequence = reverse_sequence(sequence)[-20:]  # Retain only the latter 20 characters
        reversed_pam = process_pam(pam)
        sequence_info[1] = reversed_sequence
        sequence_info[2] = reversed_pam

    # Reconstruct the modified fourth column
    modified_column_4 = ','.join(sequence_info)
    columns[3] = modified_column_4
    return "\t".join(columns)

def main():
    parser = argparse.ArgumentParser(description='Process sgRNA BED file.')
    parser.add_argument('input_file', type=str, help='Path to the input BED file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    parser.add_argument('--reverse', action='store_false', default=True, help='Reverse the sequence in the fourth column (default: True)')
    
    args = parser.parse_args()
    
    with open(args.input_file, 'r') as infile, open(args.output_file, 'w') as outfile:
        header = next(infile)  # Read the header line
        outfile.write(header)  # Write the header to the output file
        for line in infile:
            modified_line = process_line(line, args.reverse)
            if modified_line:  # Only write if the line was correctly processed
                outfile.write(f"{modified_line}\n")
    
    print("Processing complete. The output has been saved to", args.output_file)

if __name__ == "__main__":
    main()