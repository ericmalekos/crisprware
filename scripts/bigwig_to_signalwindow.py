#!/usr/bin/env python3

'''
Usage:
./biwig_to_signalwindow.py <reference.bed> <signal.bw> <output.bed> -w 150

'''
import argparse
import pyBigWig

# Set up argument parsing
parser = argparse.ArgumentParser(description="Find a specified window with the highest average signal in a BigWig file for each entry in a BED file.")
parser.add_argument("bed_file", help="Path to the input BED file")
parser.add_argument("bigwig_file", help="Path to the input BigWig file")
parser.add_argument("output_bed_file", help="Path to the output BED file")
parser.add_argument("-w", "--window_size", type=int, default=200, help="Size of the window to analyze (default: 100)")

# Parse arguments
args = parser.parse_args()

# Function to find the window with the highest average signal
def find_highest_signal_window(bw, chrom, start, end, window_size):
    max_avg = float('-inf')
    max_start = start
    for i in range(start, end - window_size + 1):
        window_end = i + window_size
        avg_signal = bw.stats(chrom, i, window_end, type='mean')[0]
        if avg_signal > max_avg:
            max_avg = avg_signal
            max_start = i
    return max_start, max_start + window_size, max_avg

# Open the BigWig file
bw = pyBigWig.open(args.bigwig_file)

with open(args.bed_file, 'r') as bed, open(args.output_bed_file, 'w') as out_bed:
    for line in bed:
        chrom, start, end, *_ = line.strip().split()
        start, end = int(start), int(end)
        
        # Find the highest signal window within the BED entry
        highest_start, highest_end, _ = find_highest_signal_window(bw, chrom, start, end, args.window_size)
        
        # Write the highest signal window to the output BED file
        out_bed.write(f"{chrom}\t{highest_start}\t{highest_end}\n")

# Close the BigWig file
bw.close()

print(f"Process completed. The output is saved in '{args.output_bed_file}'. Window size used: {args.window_size} bp.")