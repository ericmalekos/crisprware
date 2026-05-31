#!/usr/bin/env python3

"""
This is a wrapper script for building off-target indices.

Two engines are supported and may be built together:
  * crispr-ots  - fast exact-mismatch scanner (default)
  * guidescan2  - supports RNA/DNA bulges

index_genome -f ./test_data/test.fasta
index_genome -f ./test_data/test.fasta --indexer crispr-ots guidescan2
"""

from typing import List, Optional

import argparse
import subprocess
from crisprware.utils.gtf_bed_processing_functions import merge_targets
from crisprware.utils.dna_sequence_functions import subset_fasta_with_bed
from crisprware.utils.utility_functions import create_output, decompress_gzip_if_needed, remove_file

# Each engine ships its own binary but shares the same `index --index <out> <fasta>`
# CLI. The extension keeps each engine's index in its own `<name>_<ext>/` directory
# so both can be built for the same FASTA without their files colliding, and so the
# downstream specificity column (`specificity_<name>_<ext>`) records which engine
# produced it.
INDEXER_EXECUTABLE = {"crispr-ots": "crispr-ots", "guidescan2": "guidescan"}
INDEXER_EXTENSION = {"crispr-ots": "crisprots", "guidescan2": "guidescan2"}


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add index_genome arguments to the given parser."""
    parser.add_argument(
        "-f", "--fasta", type=str, help="FASTA file to use as a reference for index creation.", required=True
    )

    parser.add_argument(
        "--indexer",
        type=str,
        nargs="+",
        choices=["crispr-ots", "guidescan2"],
        default=["crispr-ots"],
        help="Off-target index engine(s) to build. 'crispr-ots' (fast exact-mismatch \
            scanner; default) and/or 'guidescan2' (slower, but supports RNA/DNA bulges). \
            Pass both (e.g. '--indexer crispr-ots guidescan2') to build one index of each \
            type; score_guides can then score against either or both. Each engine writes \
            to its own '<name>_<engine>' directory. [default: crispr-ots]",
    )

    parser.add_argument(
        "-k",
        "--locations_to_keep",
        help="List of BED/GTF files with coordinates to use \
        for index creation. These locations will be used for \
        off-target scoring. If multiple files are passed, coordinates \
        will be merged with a union operation. \
        Leave empty to use entire fasta.",
        type=str,
        default="",
        nargs="*",
    )

    parser.add_argument(
        "--feature",
        type=str,
        help="For any GTF/GFF in '--locations_to_keep', only this \
            feature will be used for determining appropriate sgRNA. \
            The feature should match an entry in the third column of \
            the GTF/GFF. [default: 'transcript']",
        default="transcript",
    )

    parser.add_argument(
        "-w",
        "--context_window",
        nargs=2,
        type=int,
        help="Pass two, space-separated, integers to specifiy the \
            nucleotide window around the --locations_to_keep \
            '<upstream>' '<downstream>'. This can be used to \
            expand the window around the final intervals e.g. \
            '-w 1000 1500' expands chr1 2000 3500 -> chr1 1000 5000\n \
            Good for CRISPRi/a \
            [default: 20 20]",
        default=[20, 20],
    )

    parser.add_argument(
        "-o", "--output_directory", help="Path to output. [default: current directory]", type=str, default=""
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Builds crispr-ots and/or Guidescan2 off-target indices.")
    add_arguments(parser)
    return parser.parse_args()


def guideScanIndex(fasta: str, output: str, indexer: str = "crispr-ots") -> int:
    """
    Build an off-target index for a FASTA file using the selected engine.

    Parameters:
        fasta (str): The path to the input fasta file.
        output (str): The path/prefix where the index will be written.
        indexer (str): Which engine to use, 'crispr-ots' (default) or 'guidescan2'.
            Both expose the same 'index --index <out> <fasta>' CLI; only the binary
            name differs.

    Returns:
        int: A success code (1) indicating that the index was generated successfully.
    """

    cmd = [INDEXER_EXECUTABLE[indexer], "index", "--index", output, fasta]

    subprocess.run(cmd, check=True)

    return 1


def main(args: Optional[argparse.Namespace] = None) -> None:

    if args is None:
        args = parse_arguments()

    # Decompress the file if it's gzipped
    # original_fasta = args.fasta
    args.fasta, was_gzipped = decompress_gzip_if_needed(args.fasta)

    # The merged-interval bed and subset fasta are engine-independent, so build them
    # once (if requested) and index the result with every selected engine.
    build_fasta = args.fasta

    if args.locations_to_keep:
        prep_path, _ = create_output(args.fasta, outdir=args.output_directory, extension="subset", tmp=False)

        bed_output_path = prep_path + "_merged.bed"
        locs_to_keep = merge_targets(
            args.locations_to_keep, gtf_feature=args.feature, operation="merge", window=args.context_window
        )
        print(f"\n\tSaving merged interval bed to {bed_output_path}")
        locs_to_keep.saveas(bed_output_path)

        fasta_output_path = prep_path + ".fasta"
        print(f"\n\tSaving subset fasta to {fasta_output_path}")
        subset_fasta_with_bed(fasta_path=args.fasta, bed_path=bed_output_path, output_fasta_path=fasta_output_path)

        build_fasta = fasta_output_path

    for indexer in args.indexer:
        index_output_path, _ = create_output(
            args.fasta, outdir=args.output_directory, extension=INDEXER_EXTENSION[indexer]
        )
        print(f"\n\tBuilding {indexer} index from {build_fasta}")
        print(f"\tSaving index to {index_output_path}")
        guideScanIndex(build_fasta, index_output_path, indexer)

    remove_file(f"{build_fasta}.forward.dna")
    remove_file(f"{build_fasta}.reverse.dna")

    if was_gzipped:
        remove_file(args.fasta)


if __name__ == "__main__":
    main()
