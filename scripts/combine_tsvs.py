#!/usr/bin/env python3

'''
Usage:
combine_tsvs.py file1.tsv file2.tsv file3.tsv --output output.tsv


'''
import pandas as pd
import argparse
import sys

def concatenate_tsvs(tsv_files):
    # Read the first TSV to get the column headers
    try:
        base_df = pd.read_csv(tsv_files[0], sep='\t')
    except Exception as e:
        print(f"Error reading {tsv_files[0]}: {e}")
        sys.exit(1)

    # Column headers from the first file
    base_columns = base_df.columns

    # List to hold all DataFrames including the first one
    dfs = [base_df]

    # Iterate over the remaining TSV files
    for tsv_file in tsv_files[1:]:
        try:
            # Read current TSV file
            df = pd.read_csv(tsv_file, sep='\t')

            # Check if the column headers match the first file
            if not df.columns.equals(base_columns):
                print(f"Column headers do not match in {tsv_file}. Exiting.")
                sys.exit(1)

            # If headers match, append the DataFrame to the list
            dfs.append(df)
        except Exception as e:
            print(f"Error reading {tsv_file}: {e}")
            sys.exit(1)

    # Concatenate all DataFrames in the list
    concatenated_df = pd.concat(dfs, ignore_index=True)

    # Return the concatenated DataFrame
    return concatenated_df

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Concatenate TSV files with matching headers.")
    parser.add_argument("tsv_files", nargs='+', help="Space-separated list of TSV files to concatenate.")
    parser.add_argument("--output", required=True, help="Output file name for the concatenated TSV.")

    # Parse command line arguments
    args = parser.parse_args()

    # Concatenate TSVs
    concatenated_df = concatenate_tsvs(args.tsv_files)

    # Save the concatenated DataFrame to a new TSV file
    concatenated_df.to_csv(args.output, sep='\t', index=False)
    print(f"Concatenated TSV saved to {args.output}")

if __name__ == "__main__":
    main()