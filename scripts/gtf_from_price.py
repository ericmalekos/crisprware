#!/usr/bin/env python3

import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter price TSV files and convert them to GTF format.")
    parser.add_argument("-i", "--input_tsv", required=True, help="Path to the input price TSV file")
    parser.add_argument("-g", "--input_gtf", required=True, help="Path to the input GTF file to be used as a reference")
    parser.add_argument("-o", "--output_gtf", required=True, help="Path to output the new GTF file")
    parser.add_argument("-p", "--min_p_value", type=float, default=1.0, help="Minimum p value for filtering")
    parser.add_argument("--min_aalen", type=int, default=1, help="Minimum amino acid length")
    parser.add_argument("--tis_type", help="Tis Type to filter")
    parser.add_argument("--start_codon", help="start codon to filter")

    return parser.parse_args()

def load_price_data(input_file):
    return pd.read_csv(input_file, sep='\t')

def load_gtf_data(gtf_file):
    gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    return pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=gtf_columns, quoting=3, dtype={'start': int, 'end': int})

def calculate_aalen(location):
    # Strip out the part before and including the colon
    clean_location = location.split(':')[-1] if ':' in location else location
    total_length = 0
    
    blocks = clean_location.split('|')
    
    for block in blocks:
        try:
            start, end = block.split('-')
            start = int(start)
            end = int(end)
            total_length += end - start + 1
        except (ValueError, IndexError) as e:
            print(f"Warning: Skipping invalid location format: '{block}'. Error: {e}")
    
    # Calculate amino acid length (nucleotides divided by 3, subtract 1 for stop codon)
    return (total_length // 3) - 1

def filter_data(df, args):
    df['AALen'] = df['Location'].apply(calculate_aalen)
    filtered_df = df[(df['p value'] <= args.min_p_value) & (df['AALen'] >= args.min_aalen)]
    
    if args.tis_type:
        filtered_df = filtered_df[filtered_df['Type'] == args.tis_type]
    if args.start_codon:
        filtered_df = filtered_df[filtered_df['Codon'] == args.start_codon]
    
    return filtered_df

def convert_to_gtf(price_df, gtf_df, output_gtf_path):
    # Extract transcript_id from the attribute column
    gtf_df['transcript_id'] = gtf_df['attribute'].str.extract('transcript_id "([^"]+)"')
    gtf_df['gene_id'] = gtf_df['attribute'].str.extract('gene_id "([^"]+)"')

    # Create a mapping from GTF transcript_id to price_df Id
    transcript_id_mapping = {row['Id'].split('_')[0]: row['Id'] for _, row in price_df.iterrows()}

    with open(output_gtf_path, 'w') as out_file:
        for _, price_row in price_df.iterrows():
            chrom_strand, clean_location = price_row['Location'].split(':', 1)
            chrom = 'chr' + chrom_strand[:-1] if not chrom_strand[:-1].startswith("chr") else chrom_strand[:-1]
            strand = chrom_strand[-1]
            tscript_id = price_row['Id']
            gene_id = price_row['Gene']

            # Find corresponding transcript in GTF
            transcript_gtf = gtf_df[(gtf_df['transcript_id'] == tscript_id.split('_')[0]) & (gtf_df['feature'] == 'transcript')]
            
            if not transcript_gtf.empty:
                # Write transcript entry
                transcript_entry = transcript_gtf.iloc[0].copy()
                transcript_entry['attribute'] = f'gene_id "{gene_id}"; transcript_id "{tscript_id}";'
                out_file.write('\t'.join(map(str, transcript_entry.iloc[:9])) + '\n')

                # Find and write exon entries
                exon_entries = gtf_df[(gtf_df['transcript_id'] == tscript_id.split('_')[0]) & (gtf_df['feature'] == 'exon')]
                for _, exon in exon_entries.iterrows():
                    exon = exon.copy()
                    exon['attribute'] = f'gene_id "{gene_id}"; transcript_id "{tscript_id}";'
                    out_file.write('\t'.join(map(str, exon.iloc[:9])) + '\n')

            # Write CDS entries
            for block in clean_location.split('|'):
                start, end = map(int, block.split('-'))
                start += 1  # GTF is 1-based
                attributes = f'gene_id "{gene_id}"; transcript_id "{tscript_id}";'
                gtf_line = f'{chrom}\tPrice\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attributes}\n'
                out_file.write(gtf_line)

def main():
    args = parse_arguments()
    price_df = load_price_data(args.input_tsv)
    gtf_df = load_gtf_data(args.input_gtf)
    filtered_df = filter_data(price_df, args)
    convert_to_gtf(filtered_df, gtf_df, args.output_gtf)

if __name__ == "__main__":
    main()
