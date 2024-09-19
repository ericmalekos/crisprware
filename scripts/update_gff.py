#!/usr/bin/env python3

import sys
import argparse

def process_gff(input_file, output_file):
    gene_id_map = {}
    transcript_id_map = {}

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue

            attributes = dict(item.split('=') for item in fields[8].split(';'))

            if fields[2] == 'gene':
                gene_id = attributes['ID']
                gene_id_map[gene_id] = gene_id
                attributes['gene_id'] = gene_id
            elif fields[2] == 'mRNA':
                fields[2] = 'transcript'
                transcript_id = attributes['ID']
                parent_gene = attributes['Parent']
                gene_id = gene_id_map.get(parent_gene, parent_gene)
                transcript_id_map[transcript_id] = (gene_id, transcript_id)
                attributes['gene_id'] = gene_id
                attributes['transcript_id'] = transcript_id
            else:
                parent = attributes.get('Parent', '')
                if parent in transcript_id_map:
                    gene_id, transcript_id = transcript_id_map[parent]
                    attributes['gene_id'] = gene_id
                    attributes['transcript_id'] = transcript_id

            # Reconstruct the attributes string
            new_attributes = ';'.join(f"{k}={v}" for k, v in attributes.items())
            fields[8] = new_attributes

            outfile.write('\t'.join(fields) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Process GFF file to match desired format.')
    parser.add_argument('input', help='Input GFF file')
    parser.add_argument('output', help='Output GFF file')
    args = parser.parse_args()

    process_gff(args.input, args.output)

if __name__ == '__main__':
    main()