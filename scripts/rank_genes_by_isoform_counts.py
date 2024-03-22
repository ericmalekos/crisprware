#!/usr/bin/env python3
'''
    ./rank_genes_by_isoform_counts.py ../../AnnotationFiles/gencode.vM31.primary_assembly.annotation.gtf --protein_coding > test.tsv

'''

import argparse
import csv

def count_isoforms(gtf_file, protein_coding_only):
    """
    Counts the number of unique isoforms for each gene in the provided GTF file,
    with an option to include only protein-coding transcripts. Also captures the gene name.
    """
    # Initialize a dictionary to hold the count of isoforms and gene names for each gene.
    gene_info = {}

    with open(gtf_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if row[0].startswith('#') or len(row) < 9:
                # Skip header lines and incomplete rows
                continue

            attributes = row[8]
            # Parse the attributes column to extract gene_id, transcript_id, gene_type, and gene_name
            attr_dict = {attr.strip().split(' ')[0]: attr.strip().split(' ')[1].strip('"') for attr in attributes.split(';') if attr}

            gene_id = attr_dict.get('gene_id', None)
            transcript_id = attr_dict.get('transcript_id', None)
            gene_name = attr_dict.get('gene_name', None)
            gene_type = attr_dict.get('gene_type', None) if protein_coding_only else 'protein_coding'

            if gene_id and transcript_id and gene_name and gene_type == 'protein_coding':
                if gene_id not in gene_info:
                    gene_info[gene_id] = {'isoforms': set(), 'gene_name': gene_name}
                # Add the transcript_id to the set of isoforms for this gene
                gene_info[gene_id]['isoforms'].add(transcript_id)

    # Convert the sets to counts and include gene names
    gene_isoform_count = {gene: {'count': len(info['isoforms']), 'gene_name': info['gene_name']} for gene, info in gene_info.items()}

    return gene_isoform_count

def main():
    parser = argparse.ArgumentParser(description="Rank genes by the number of unique isoforms from a GTF file, with an option to include only protein-coding transcripts. Outputs gene ID, gene name, and the number of isoforms.")
    parser.add_argument("gtf_file", help="Path to the GTF file")
    parser.add_argument("--protein_coding", action="store_true", help="Include only protein-coding transcripts in the count")
    
    args = parser.parse_args()

    isoform_counts = count_isoforms(args.gtf_file, args.protein_coding)

    # Rank genes by the number of unique isoforms
    ranked_genes = sorted(isoform_counts.items(), key=lambda x: x[1]['count'], reverse=True)

    # Print the ranked list
    print("GeneID\tGeneName\tNumber of Isoforms")
    for gene, info in ranked_genes:
        print(f"{gene}\t{info['gene_name']}\t{info['count']}")

if __name__ == "__main__":
    main()