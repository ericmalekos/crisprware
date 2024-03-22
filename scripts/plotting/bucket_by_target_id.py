#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
from itertools import accumulate
import matplotlib.cm as cm
import numpy as np

'''
        ./bucket_by_target_id.py ../ChIP_Methods/RankedSgRNAs/RankedSgRNA.tsv \
        specificity_encode_metagene specificity_Hg38_index \
        --title "CEBPB K562 Targets" --output CEBPB_K562.png
'''


def count_targets_by_specificity(tsv_file, column_names):
    counts = {col: defaultdict(int) for col in column_names}
    with open(tsv_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            target_id = row['target_id']
            for column_name in column_names:
                try:
                    value = float(row[column_name])
                    if value > 0.2:
                        counts[column_name][target_id] += 1
                except ValueError:
                    # Skip rows where the column does not have a float value
                    continue

    return counts

def bucket_counts(counts):
    bucket_labels = ['0', '1', '2', '3-4', '5-9', '10-19', '20-49', '>=50']
    bucket_ranges = [0, 1, 2, range(3, 5), range(5, 10), range(10, 20), range(20, 50), float('inf')]
    buckets = {col: {label: 0 for label in bucket_labels} for col in counts.keys()}

    for col, col_counts in counts.items():
        for count in col_counts.values():
            for label, rng in zip(bucket_labels, bucket_ranges):
                if (isinstance(rng, range) and count in rng) or count == rng:
                    buckets[col][label] += 1
                    break
                elif rng == float('inf') and count >= 50:
                    buckets[col]['>=50'] += 1
                    break

    print(buckets)
    return buckets

def plot_buckets(all_buckets, output_file, plot_title):
    column_names = list(all_buckets.keys())
    bucket_labels = list(next(iter(all_buckets.values())).keys())
    
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define the colors for the buckets
    color_map = cm.get_cmap('Reds', len(bucket_labels))

    # Plotting each column as a stacked bar
    width = 0.65  # the width of the bars
    column_indices = np.arange(len(column_names))  # the column locations on the x-axis

    # Initialize the bottom for the stacked bars for each column
    bottoms = np.zeros(len(column_names))

    for i, bucket_label in enumerate(bucket_labels):
        # Get the values for this bucket across all columns
        bucket_values = np.array([all_buckets[column_name][bucket_label] for column_name in column_names])

        # Normalize the bucket values to sum to 1 within each column
        column_totals = np.array([sum(all_buckets[column_name].values()) for column_name in column_names])
        normalized_values = bucket_values / column_totals  # Normalize the values for each column

        # Get the color for this bucket
        color = color_map(i / len(bucket_labels))

        # Create the stacked bar for this bucket across all columns
        ax.bar(column_indices, normalized_values, width, bottom=bottoms, color=color, label=bucket_label)

        # Update bottoms for the next bucket
        bottoms += normalized_values

    ax.set_ylabel('Fraction of motifs')
    ax.set_title(plot_title)
    ax.set_xticks(column_indices)
    ax.set_xticklabels(column_names)
    ax.legend(title='Bucket Ranges', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')
    #plt.show()

def main():
    parser = argparse.ArgumentParser(description="Count occurrences of target_id entries based on user-defined specificity columns, plot the results, and save the figure.")
    parser.add_argument("tsv_file", help="Path to the TSV file")
    parser.add_argument("column_names", nargs='+', help="Space-separated list of column names to filter by")
    parser.add_argument("--title", default="GuideScan score distribution", help="Title for the plot")
    parser.add_argument("--output", default="plot.png", help="Output file name for the plot")
    args = parser.parse_args()

    counts = count_targets_by_specificity(args.tsv_file, args.column_names)
    buckets = bucket_counts(counts)
    
    plot_buckets(buckets, args.output, args.title)

if __name__ == "__main__":
    main()