#!/usr/bin/env python3

'''
Usage:
./bigwig_to_signalwindow.py <reference.bed> <signal.bw> <output.bed> -w 150

'''
import argparse
import pyBigWig

# Set up argument parsing
parser = argparse.ArgumentParser(description="Find a specified window with the highest average signal in a BigWig file for each entry in a BED file.")
parser.add_argument("bed_file", help="Path to the input BED file")
parser.add_argument("bigwig_file", help="Path to the input BigWig file")
parser.add_argument("output_bed_file", help="Path to the output BED file")
parser.add_argument("-w", "--window_size", type=int, default=200, help="Size of the window to analyze (default: 200)")

# Parse arguments
args = parser.parse_args()

# Function to find the window with the highest average signal
def find_highest_signal_window(bw, chrom, start, end, window_size, chrom_size):
    start = max(0, start)
    end = min(chrom_size, end)
    
    if end - start < window_size:
        return None, None, None

    max_avg = float('-inf')
    max_start = start
    valid_window_found = False

    for i in range(start, end - window_size + 1):
        window_end = i + window_size
        try:
            avg_signal = bw.stats(chrom, i, window_end, type='mean')[0]
            if avg_signal is not None:
                valid_window_found = True
                if avg_signal > max_avg:
                    max_avg = avg_signal
                    max_start = i
        except RuntimeError:
            continue  # Skip this window if there's an error

    if not valid_window_found:
        return None, None, None
    return max_start, max_start + window_size, max_avg

# Open the BigWig file
bw = pyBigWig.open(args.bigwig_file)

# Get the list of chromosomes and their sizes from the BigWig file
bigwig_chroms = bw.chroms()

print("Processing BED file: {}".format(args.bed_file))
print("Using BigWig file: {}".format(args.bigwig_file))
print("Output will be saved to: {}".format(args.output_bed_file))
print("Window size: {} bp".format(args.window_size))
print("Number of chromosomes in BigWig file: {}".format(len(bigwig_chroms)))

processed_entries = 0
skipped_entries = 0

with open(args.bed_file, 'r') as bed, open(args.output_bed_file, 'w') as out_bed:
    for line in bed:
        chrom, start, end, *rest = line.strip().split()
        start, end = int(start), int(end)
        
        # Check if the chromosome is in the BigWig file
        if chrom not in bigwig_chroms:
            skipped_entries += 1
            continue
        
        # Find the highest signal window within the BED entry
        highest_start, highest_end, max_avg = find_highest_signal_window(bw, chrom, start, end, args.window_size, bigwig_chroms[chrom])
        
        # Write the highest signal window to the output BED file
        if highest_start is not None:
            out_bed.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, highest_start, highest_end, max_avg, "\t".join(rest)))
        else:
            out_bed.write("{}\t{}\t{}\tNo_valid_data\t{}\n".format(chrom, start, end, "\t".join(rest)))
        
        processed_entries += 1

# Close the BigWig file
bw.close()

print("Process completed.")
print("Processed entries: {}".format(processed_entries))
print("Skipped entries (chromosome not in BigWig): {}".format(skipped_entries))
print("Output saved in: {}".format(args.output_bed_file))