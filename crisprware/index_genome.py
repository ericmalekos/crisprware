#!/usr/bin/env python3

"""
Wrapper for building a crispr-ots off-target index. The enzyme is defined by a
PAM motif, protospacer length, and PAM orientation (5'/3') — no enzyme presets.

  index_genome -f ./test_data/test.fasta --pam NGG -l 20            # SpCas9
  index_genome -f ./genome.fa --pam TTTV -l 23 --pam_5_prime --bin_width 14  # Cas12a
"""

from typing import Optional

import argparse
import subprocess
from crisprware.utils.gtf_bed_processing_functions import merge_targets
from crisprware.utils.dna_sequence_functions import subset_fasta_with_bed
from crisprware.utils.utility_functions import create_output, decompress_gzip_if_needed, remove_file


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add index_genome arguments to the given parser."""
    parser.add_argument(
        "-f", "--fasta", type=str, help="FASTA file to use as a reference for index creation.", required=True
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

    # --- Enzyme definition for the crispr-ots index (PAM/length/orientation,
    #     no presets) ---
    parser.add_argument(
        "-p",
        "--pam",
        type=str,
        required=True,
        help="PAM motif as an IUPAC string, e.g. NGG (SpCas9), NAG, TTTV or TTTN \
            (Cas12a). With --protospacer_length and --pam_5_prime this fully \
            defines the indexed off-target sites.",
    )

    parser.add_argument(
        "-l",
        "--protospacer_length",
        type=int,
        required=True,
        help="Protospacer (spacer) length in bases, e.g. 20 for SpCas9 or 23 for Cas12a.",
    )

    parser.add_argument(
        "--pam_5_prime",
        action="store_true",
        default=False,
        help="Set if the PAM is 5' of the protospacer (Cas12a-style). Omit for a \
            3' PAM (SpCas9-style). [default: False]",
    )

    parser.add_argument(
        "--bin_width",
        type=int,
        default=None,
        help="crispr-ots bin-prefix width (1-15); higher = faster GPU scan but a \
            larger index. Sweet spots: ~9-10 SpCas9 full-human, 14 for Cas12a \
            full-human on a 24 GB GPU. [default: engine default]",
    )

    parser.add_argument(
        "-o", "--output_directory", help="Path to output. [default: current directory]", type=str, default=""
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Builds a crispr-ots off-target index for a given PAM / protospacer length / orientation."
    )
    add_arguments(parser)
    return parser.parse_args()


def guideScanIndex(
    fasta: str,
    output: str,
    pam: str,
    protospacer_length: int,
    pam_5_prime: bool = False,
    bin_width: Optional[int] = None,
) -> int:
    """
    Build a crispr-ots off-target index for a FASTA, for the enzyme defined by
    the PAM motif, protospacer length, and PAM orientation (no enzyme presets).

    Parameters:
        fasta (str): Path to the input FASTA.
        output (str): Output index prefix.
        pam (str): IUPAC PAM motif (e.g. NGG, TTTV).
        protospacer_length (int): Spacer length in bases (e.g. 20, 23).
        pam_5_prime (bool): True if the PAM is 5' of the protospacer (Cas12a).
        bin_width (int | None): crispr-ots bin-prefix width (1-15), or None for the default.

    Returns:
        int: A success code (1) indicating that the index was generated successfully.
    """

    cmd = ["crispr-ots", "index", "--index", output, "--pam", pam, "-l", str(protospacer_length)]
    if pam_5_prime:
        cmd.append("--pam-5-prime")
    if bin_width is not None:
        cmd += ["-w", str(bin_width)]
    cmd.append(fasta)

    subprocess.run(cmd, check=True)

    return 1


def main(args: Optional[argparse.Namespace] = None) -> None:

    if args is None:
        args = parse_arguments()

    # Decompress the file if it's gzipped
    # original_fasta = args.fasta
    args.fasta, was_gzipped = decompress_gzip_if_needed(args.fasta)

    index_output_path, _ = create_output(args.fasta, outdir=args.output_directory, extension="gscan2")

    if not args.locations_to_keep:
        guideScanIndex(
            args.fasta, index_output_path, args.pam, args.protospacer_length, args.pam_5_prime, args.bin_width
        )
    else:
        bed_output_path, _ = create_output(args.fasta, extension="gscan2", tmp=False)
        bed_output_path += "_merged.bed"

        locs_to_keep = merge_targets(
            args.locations_to_keep, gtf_feature=args.feature, operation="merge", window=args.context_window
        )

        print(f"\n\tSaving merged interval bed to {bed_output_path}")

        locs_to_keep.saveas(bed_output_path)

        fasta_output_path, _ = create_output(args.fasta, extension="gscan2", tmp=False)
        fasta_output_path += "_subset.fasta"

        print(f"\n\tSaving subset fasta to {fasta_output_path}")

        subset_fasta_with_bed(fasta_path=args.fasta, bed_path=bed_output_path, output_fasta_path=fasta_output_path)

        print("\tBuilding Index from " + fasta_output_path)
        print("\tSaving Index to " + index_output_path)
        guideScanIndex(
            fasta_output_path,
            index_output_path.strip("_"),
            args.pam,
            args.protospacer_length,
            args.pam_5_prime,
            args.bin_width,
        )

        args.fasta = fasta_output_path

    remove_file(f"{args.fasta}.forward.dna")
    remove_file(f"{args.fasta}.reverse.dna")

    if was_gzipped:
        remove_file(args.fasta)


if __name__ == "__main__":
    main()
