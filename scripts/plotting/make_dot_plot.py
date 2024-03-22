#!/usr/bin/env python3

'''
    ./make_dot_plot.py ../../ChIP_Methods/RankedSgRNAs/RankedSgRNA.tsv \
    specificity_encode_metagene specificity_Hg38_index --sample_size 1000 \
    --title "CEBPB_sgRNA_targets" --output CEBPB_dotplot.png

'''
import argparse
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

def count_targets_by_specificity(tsv_file, column_names):
    target_counts = {col: defaultdict(int) for col in column_names}
    with open(tsv_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            target_id = row['target_id']
            if any(row[col] == "-1" for col in column_names):
                continue
            for column_name in column_names:
                try:
                    value = float(row[column_name])
                    if value > 0.2:
                        target_counts[column_name][target_id] += 1
                except ValueError:
                    continue
    return target_counts

def perform_mann_whitney_u_test(values1, values2):
    stat, p = mannwhitneyu(values1, values2)
    return p

def make_dot_and_violin_plot(target_counts, column_names, sample_size, output_file, plot_title):
    full_data_values = [list(target_counts[col].values()) for col in column_names]
    p_value = perform_mann_whitney_u_test(full_data_values[0], full_data_values[1])

    data_for_plotting = []
    for values in full_data_values:
        if sample_size and sample_size < len(values):
            downsampled_values = np.random.choice(values, size=sample_size, replace=False)
        else:
            downsampled_values = values
        data_for_plotting.append(downsampled_values)
    
    fig, ax = plt.subplots()
    parts = ax.violinplot(data_for_plotting, showmeans=False, showmedians=False, showextrema=False)
    medians = [np.median(values) for values in data_for_plotting]
    for i, values in enumerate(data_for_plotting):
        x_values = np.random.normal(i + 1, 0.04, size=len(values))
        ax.scatter(x_values, values, alpha=0.6, s=10)
        ax.plot([i + 0.75, i + 1.25], [medians[i]]*2, color='k', linestyle='-', linewidth=2)
    
    ax.annotate(f'p={p_value:.2e}', xy=(0.5, -0.15), xycoords='axes fraction', ha='center', va='top')
    
    ax.set_ylabel('Count of values > 0.2')
    ax.set_title(plot_title)
    ax.set_xticks(np.arange(1, len(column_names) + 1))
    ax.set_xticklabels(column_names)

    if output_file:
        plt.savefig(output_file, bbox_inches='tight')
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Create a dot and violin plot for specified columns based on value counts greater than 0.2, excluding rows with -1, perform Mann-Whitney U test on the full dataset, and save the plot to a file.")
    parser.add_argument("tsv_file", help="Path to the TSV file")
    parser.add_argument("column_names", nargs='+', help="Column names for which to generate the plot")
    parser.add_argument("--sample_size", type=int, help="Maximum number of samples to include in the plot for each column (optional)")
    parser.add_argument("--output", type=str, help="The file path to save the plot image")
    parser.add_argument("--title", type=str, default="Dot and Violin Plot of Counts by Column", help="Title for the plot")

    args = parser.parse_args()

    target_counts = count_targets_by_specificity(args.tsv_file, args.column_names)
    make_dot_and_violin_plot(target_counts, args.column_names, args.sample_size, args.output, args.title)

if __name__ == "__main__":
    main()