#!/usr/bin/env python3

'''
    This is a wrapper script for generating Guidescan2 indices

    index_genome -f ./test_data/test.fasta
'''
import argparse
import subprocess
from utils.gtf_bed_processing_functions import merge_targets
from utils.dna_sequence_functions import subset_fasta_with_bed
from utils.utility_functions import create_output, decompress_gzip_if_needed, remove_file


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Builds Guidescan2 Index."
    )

    parser.add_argument(
        "-f", "--fasta",
        type=str,
        help="FASTA file to use as a reference for index creation.",
        required=True
    )

    parser.add_argument(
        "-k", "--locations_to_keep",
        help="List of BED/GTF files with coordinates to use \
        for index creation. These locations will be used for \
        off-target scoring. If multiple files are passed, coordinates \
        will be merged with a union operation. \
        Leave empty to use entire fasta.",
        type=str,
        default="",
        nargs='*'
    )

    parser.add_argument(
        "--feature",
        type=str,
        help="For any GTF/GFF in '--locations_to_keep', only this \
            feature will be used for determining appropriate sgRNA. \
            The feature should match an entry in the third column of \
            the GTF/GFF. [default: 'transcript']",
        default="transcript"
    )

    parser.add_argument(
        "-w", "--context_window",
        nargs=2,
        type=int,
        help="Pass two, space-separated, integers to specifiy the \
            nucleotide window around the --locations_to_keep \
            '<upstream>' '<downstream>'. This can be used to \
            expand the window around the final intervals e.g. \
            '-w 1000 1500' expands chr1 2000 3500 -> chr1 1000 5000\n \
            Good for CRISPRi/a \
            [default: 20 20]",
        default=[20,20]
    )

    parser.add_argument(
        "-o", "--output_directory",
        help="Path to output. [default: current directory]",
        type=str,
        default=""
    )

    return parser.parse_args()

def guideScanIndex(fasta, output):
    """
    This function generates a Guidescan index for a given fasta file.

    Parameters:
        fasta (str): The path to the input fasta file.
        output (str): The path where the output index will be saved.

    Returns:
        int: A success code (1) indicating that the index was generated successfully.
    """

    cmd = [
        'guidescan',
        'index',
        '--index',
        output,
        fasta
    ]

    subprocess.run(cmd, check=True)

    return 1



def main():

    args = parse_arguments()

    # Decompress the file if it's gzipped
    #original_fasta = args.fasta
    args.fasta, was_gzipped = decompress_gzip_if_needed(args.fasta)


    index_output_path, _ = create_output(args.fasta, outdir=args.output_directory, extension="gscan2")

    if not args.locations_to_keep:
        guideScanIndex(args.fasta, index_output_path)
    else:
    
        bed_output_path, _ = create_output(args.fasta, extension="gscan2", tmp=False)
        bed_output_path += "_merged.bed"

        locs_to_keep = merge_targets(args.locations_to_keep, gtf_feature=args.feature,
                                     operation="merge", window=args.context_window)

        print(f"\n\tSaving merged interval bed to {bed_output_path}")

        locs_to_keep.saveas(bed_output_path)    

        fasta_output_path, _ = create_output(args.fasta, extension="gscan2", tmp=False)
        fasta_output_path += "_subset.fasta"

        print(f"\n\tSaving subset fasta to {fasta_output_path}")

        subset_fasta_with_bed(fasta_path=args.fasta,
                              bed_path=bed_output_path,
                              output_fasta_path=fasta_output_path)

        print("\tBuilding Index from " + fasta_output_path)
        print("\tSaving Index to " + index_output_path)
        guideScanIndex(fasta_output_path, index_output_path.strip("_"))

        args.fasta = fasta_output_path
    
    remove_file(f"{args.fasta}.forward.dna")
    remove_file(f"{args.fasta}.reverse.dna")

    if was_gzipped:
        remove_file(args.fasta)

if __name__ == "__main__":
    main()
