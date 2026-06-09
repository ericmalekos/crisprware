#!/usr/bin/env python3

"""
This modifies "generate_kmers.py" script from the Guidescan repo

EX Usage:
    generate_guides --pam NGG  --sgRNA_length 20  -f ./chr21Index/chr21.fa --discard_poly_T --threads 8

    generate_guides --sgRNA_length 20  -f ./chr21Index/chr21.fa --discard_poly_T --threads 8 --locations_to_keep ./AnnotationFiles/NCAM2.gtf -o "NCAM_"


"""

from typing import Dict, Generator, List, Optional, Tuple, Union

import argparse
import os
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor

import pysam
from pybedtools import BedTool
from crisprware.utils.dna_sequence_functions import (
    NTS,
    map_ambiguous_sequence,
    revcom,
    merge_targets,
    include_sgRNA,
    get_chromosome_boundaries,
)
from crisprware.utils.utility_functions import create_output, decompress_gzip_if_needed, remove_file


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add generate_guides arguments to the given parser."""
    parser.add_argument(
        "-f", "--fasta", type=str, help="FASTA file to use as a reference for sgRNA generation.", required=True
    )

    parser.add_argument(
        "-p",
        "--pam",
        type=str,
        help="Protospacer adjacent motif to match.\
             All IUPAC ambiguity codes are accepted as well as standard ATCG. [default: NGG]",
        default="NGG",
    )

    parser.add_argument("-l", "--sgRNA_length", help="Length of sgRNA to generate. [default: 20]", type=int, default=20)

    parser.add_argument(
        "-w",
        "--context_window",
        nargs=2,
        type=int,
        help="Pass two, space-separated, integers to specifiy the nucleotide window around the sgRNA \
            as '<upstream>' '<downstream>'. This can be used for downstream scoring, \
            For Ruleset3 use -w 4 6 to obtain an appropriate score context. [default: 4 6]",
        default=[4, 6],
    )

    parser.add_argument(
        "-5",
        "--active_site_offset_5",
        type=int,
        help="Where cut occurs relative to PAM 5' end. \
             To avoid error, use '=' sign when passing a negative number, \
                 e.g. --active_site_offset_5=-1 [default: -4]",
        default="-4",
    )

    parser.add_argument(
        "-3",
        "--active_site_offset_3",
        type=int,
        help="Where cut occurs relative to PAM 5' end. [default: -2] \
            To avoid error, use '=' sign when passing a negative number, \
            e.g. --active_site_offset_3=-3 [default: -4]",
        default="-4",
    )

    parser.add_argument(
        "-k",
        "--locations_to_keep",
        help="List of BED/GTF files with coordinates in \
            which the sgRNA desired. If the sgRNA cutsite does not intersect \
            coordinates in these files they are discarded. Leave blank \
            to keep all sgRNA. \
            e.g. atac_peak.bed genes.gtf",
        type=str,
        default="",
        nargs="*",
    )

    parser.add_argument(
        "--feature",
        type=str,
        help="For any GTF/GFF in '--locations_to_keep', only this \
            feature will be used for determining appropriate sgRNA. \
            The feature should match an entry in the third column of the GTF/GFF. \
            [default: 'exon']",
        default="exon",
    )

    parser.add_argument(
        "--join_operation",
        type=str,
        choices=["merge", "intersect"],
        help="How to treat '--locations_to_keep' if multiple files are passed. \
        Either 'merge' or 'intersect' can be used and work as described in Bedtools. \
        If 'merge',  sgRNA will be kept if its cutsite intersects an entry in ANY of the files, if \
        'intersect' the cutsite  must intersect an entry in EACH file. [default: 'intersect']",
        default="intersect",
    )

    parser.add_argument(
        "--locations_to_discard",
        help="List of BED/GTF files with coordinates where \
            sgRNA should not target. If the sgRNA cutsite intersects \
            coordinates in these files the sgRNA is discarded. Leave blank \
            to keep all sgRNA. \
            e.g. TSS.bed coding_genes.gtf",
        type=str,
        default="",
        nargs="*",
    )

    parser.add_argument("--prefix", help="Prefix to use for sgRNA identifiers. [default: None]", type=str, default="")

    parser.add_argument(
        "--gc_range",
        nargs=2,
        type=int,
        help="Pass two, space-separated, integers to specifiy the percentile \
         range of GC content e.g. '--gc_range 25 75'. [default: 0 100]",
        default=[0, 100],
    )

    parser.add_argument(
        "--discard_poly_T",
        help="Whether to discard polyT (>TTT) sgRNA. \
            Recommend True for PolIII promoters [default: False]",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "--discard_poly_G",
        help="Whether to discard polyT (>GGGG) sgRNA. [default: False]",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "--restriction_patterns",
        help="Reject sgRNA with these restriction patterns. \
            Also checks 5'flank+sgRNA+3'flank, and reverse complement, if provided. \
            For multiple values, separate by space. \
            e.g. GCGGCCGC TCTAGA CACCTGC",
        type=str,
        default="",
        nargs="*",
    )

    parser.add_argument(
        "--flank_5",
        help="include the 5' context of the lentivirus vector. \
            Used in conjunction with --restriction_patterns to \
            remove incompatible sgRNA",
        type=str,
        default="",
    )
    parser.add_argument(
        "--flank_3",
        help="include the 3' context of the lentivirus vector. \
            Used in conjunction with --restriction_patterns to \
            remove incompatible sgRNA",
        type=str,
        default="",
    )

    parser.add_argument(
        "--min_chr_length",
        help="Minimum chromosome length to consider for sgRNA generation. [default: 20]",
        type=int,
        default=20,
    )

    parser.add_argument(
        "--pam_5_prime",
        help="If the PAM is positioned 5' to the protospacer \
        set this flag, e.g. for Cas12a sgRNAs [default: False]",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "--coords_as_active_site",
        help="Whether to output bed coordinates at the active site \
        rather than the coordinates of the entire protospacer. For purposes \
        of keeping or discarding sgRNAs, overlap with the active site \
        coordinates will be used regardless [default: True]",
        default=False,
        action="store_true",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help="Number of worker processes. Parallelism scales with the number of genome "
        "windows (see --chunk_size), not the PAM count, so this can usefully exceed it. [default: 4]",
        default=4,
    )

    parser.add_argument(
        "--chunk_size",
        type=int,
        default=5_000_000,
        help="Genome window size in bp for parallel scanning; smaller = more/finer tasks "
        "(better load balance, slightly more overhead). [default: 5,000,000]",
    )

    parser.add_argument(
        "-o", "--output_directory", help="Path to output. [default: current directory]", type=str, default=""
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generates a set of sgRNA matching user specifications. "
        "Default settings are for active SpCas9 [--pam=NGG --active_site_offset_5=-4 --active_site_offset_3=-3 --sgRNA_length=20]"
    )
    add_arguments(parser)
    return parser.parse_args()


def find_sgRNA(
    args: argparse.Namespace, pam: str, chrm: str, start: int, end: int, forward: bool = True
) -> Generator[Tuple[str, int, str], None, None]:
    """
    Generator function to find sgRNA sequences within a specified genomic region.

    Parameters:
    - args (Namespace): A Namespace object containing search parameters, including:
        - context_window (tuple): A tuple (upstream, downstream) specifying the number of
          nucleotides to include in the context sequence surrounding the sgRNA.
        - sgRNA_length (int): The length of the sgRNA sequence to extract.
        - pam_5_prime (bool): If True, indicates that the PAM sequence is located at the
          5' end of the sgRNA. Otherwise, it is assumed to be at the 3' end.
    - pam (str): The PAM sequence to search for.
    - chrm (str): The chromosome (or any DNA sequence) in which to search for sgRNA sequences.
    - start (int): The start position in the chromosome for the search.
    - end (int): The end position in the chromosome for the search.
    - forward (bool, optional): Indicates the direction of the search. If True (default),
      the search is conducted in the forward direction. If False, the search is conducted
      in the reverse direction.

    Yields:
    - tuple: A tuple containing:
        - sgRNA (str): The found sgRNA sequence in uppercase.
        - position (int): The start position of the sequence, 1-based, and according to Guidescan2 expectations
        - context (str): The context sequence surrounding the sgRNA, in uppercase.

    """
    upstr, downstr = args.context_window[0], args.context_window[1]
    index = start
    k = args.sgRNA_length
    while True:
        index = chrm.find(pam, index)

        if index == -1 or index > end:
            break

        if not args.pam_5_prime:
            if forward:
                sgRNA = chrm[index - k : index]
                context = chrm[index - k - upstr : index + downstr]
                position = index - k + 1
            else:
                sgRNA = chrm[index + len(pam) : index + k + len(pam)]
                context = chrm[index - (downstr - len(pam)) : index + k + len(pam) + upstr]
                position = index + 1
        else:
            if forward:
                sgRNA = chrm[index + len(pam) : index + len(pam) + k]
                context = chrm[index - upstr : index + len(pam) + k + downstr]
                position = index + 1
            else:
                sgRNA = chrm[index - k : index]
                context = chrm[index - k - downstr : index + len(pam) + upstr]
                position = index + 1

        index += 1

        if position < 0:
            continue

        yield sgRNA.upper(), position, context.upper()


def output_bed_line(args: argparse.Namespace, chrm_name: str, sgRNA: Dict[str, Union[str, int]]) -> str:
    """
    Formats sgRNA information into a BED line for downstream analysis or visualization.

    Parameters:
    - args (Namespace): A Namespace object containing the following attributes:
        - prefix (str): A prefix string for the sgRNA identifier.
        - active_site_offset_5 (int): The offset from the sgRNA start to the 5' end of the cleavage site.
        - active_site_offset_3 (int): The offset from the sgRNA start to the 3' end of the cleavage site.
        - pam (str): The PAM sequence used in sgRNA identification.
    - chrm_name (str): The name of the chromosome where the sgRNA is located.
    - sgRNA (dict): A dictionary containing sgRNA details, including:
        - sequence (str): The nucleotide sequence of the sgRNA.
        - position (int): The start position of the sgRNA on the chromosome.
        - sense (str): The strand orientation of the sgRNA ('+' for forward and '-' for reverse).
        - length (int): The length of the sgRNA.
        - context (str): The genomic context sequence surrounding the sgRNA.

    Returns:
    - str: A string formatted as a BED line representing the sgRNA, its cleavage site, and context.
    """

    identifier = f"{args.prefix}{chrm_name}:{sgRNA['position']}:{sgRNA['sense']}"
    guideScanEntry = ",".join(
        [identifier, str(sgRNA["sequence"]), sgRNA["pam"], chrm_name, str(sgRNA["position"]), sgRNA["sense"]]
    )

    # Find the position of the Cleavage site
    if not args.pam_5_prime:
        if sgRNA["sense"] == "+":
            cleavePos1 = sgRNA["position"] + sgRNA["length"] + args.active_site_offset_5 - 1
            cleavePos2 = sgRNA["position"] + sgRNA["length"] + args.active_site_offset_3 - 1
        else:
            cleavePos1 = sgRNA["position"] - args.active_site_offset_3 + len(args.pam) - 1
            cleavePos2 = sgRNA["position"] - args.active_site_offset_5 + len(args.pam) - 1
    else:
        if sgRNA["sense"] == "+":
            cleavePos1 = sgRNA["position"] + len(args.pam) + args.active_site_offset_5 - 1
            cleavePos2 = sgRNA["position"] + len(args.pam) + args.active_site_offset_3 - 1
        else:
            cleavePos1 = sgRNA["position"] - args.active_site_offset_3 - 1
            cleavePos2 = sgRNA["position"] - args.active_site_offset_5 - 1

    bedLine = [chrm_name, str(cleavePos1), str(cleavePos2), guideScanEntry, sgRNA["context"], sgRNA["sense"]]
    return "\t".join(bedLine)


_WORKER_FASTA = None


def _init_worker(fasta_path: str) -> None:
    """ProcessPoolExecutor initializer: open the FASTA once per worker process so
    windows are read by indexed random access (pysam), never pickled."""
    global _WORKER_FASTA
    _WORKER_FASTA = pysam.FastaFile(fasta_path)


def process_window(chrom, core_start, core_end, args, pam_set, rev_pam_set, pad, tmpdir, apply_reverse):
    """Scan every PAM over one chromosome window [core_start, core_end) and write the
    matching guides to a temp file; return its path (or None if empty).

    The worker reads only its padded window slice from the shared FASTA handle (opened
    once per worker), uppercases it once, and scans all PAMs over it -- so the genome is
    never pickled and the chromosome is never re-stringified per PAM. `pad` guarantees
    the full spacer+context of any PAM anchored in the core is in-slice; find_sgRNA's
    start/end restrict scanning to the core, so each guide is emitted exactly once.
    `apply_reverse` runs the cut-site->protospacer transform here (parallel) when there
    is no region filtering; otherwise rows stay in cut-site coords for intersect/subtract.
    """
    chrom_len = _WORKER_FASTA.get_reference_length(chrom)
    slice_start = max(0, core_start - pad)
    slice_end = min(chrom_len, core_end + pad)
    seq = _WORKER_FASTA.fetch(chrom, slice_start, slice_end).upper()
    scan_start = core_start - slice_start
    scan_end = core_end - slice_start - 1  # find_sgRNA breaks at index > end
    out = []
    for pam in pam_set + rev_pam_set:
        context_length = args.sgRNA_length + args.context_window[0] + args.context_window[1]
        if args.pam_5_prime:
            context_length += len(pam)
        forward = pam in pam_set
        for sgRNA, pos, context in find_sgRNA(
            args=args, pam=pam, chrm=seq, start=scan_start, end=scan_end, forward=forward
        ):
            if len(sgRNA) != args.sgRNA_length or len(context) != context_length:
                continue
            if not all(nt in NTS for nt in context):
                continue
            position = pos + slice_start  # slice-relative -> absolute (matches the original full-chromosome position)
            if forward:
                sense = "+"
            else:
                sgRNA = revcom(sgRNA)
                context = revcom(context)
                sense = "-"
            entry = {
                "sequence": sgRNA,
                "position": position,
                "pam": pam,
                "sense": sense,
                "length": args.sgRNA_length,
                "context": context,
            }
            if include_sgRNA(args, entry):
                out_entry = dict(entry)
                out_entry["pam"] = args.pam
                line = output_bed_line(args, chrom, out_entry)
                if apply_reverse:
                    line = reverse_cut_site_offset(line, args)
                out.append(line)
    if not out:
        return None
    fd, path = tempfile.mkstemp(suffix=".bed", dir=tmpdir)
    with os.fdopen(fd, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return path


def reverse_cut_site_offset(bedline: str, args: argparse.Namespace) -> str:
    """
    Given a line in bed format, reverse the calculatations that have
    gone into calculating the cut_site
    """
    bedlist = bedline.split("\t")
    position = int(bedline.split(",")[-2])
    # print(bedlist)
    # if bedlist[-1].strip() == "+":
    #     bedlist[1] = str(int(bedlist[1]) - args.sgRNA_length - args.active_site_offset_5)
    #     bedlist[2] = str(int(bedlist[2]) - args.active_site_offset_3)
    # else:
    #     bedlist[1] = str(int(bedlist[1]) + args.active_site_offset_5)
    #     bedlist[2] = str(int(bedlist[2]) + args.sgRNA_length + args.active_site_offset_3)

    if not args.pam_5_prime:
        if bedlist[-1].strip() == "+":
            bedlist[1] = position - 1
            bedlist[2] = position + args.sgRNA_length - 1
        else:
            bedlist[1] = position + len(args.pam) - 1
            bedlist[2] = position + args.sgRNA_length + len(args.pam) - 1
    else:
        if bedlist[-1].strip() == "+":
            bedlist[1] = position + len(args.pam) - 1
            bedlist[2] = position + len(args.pam) + args.sgRNA_length - 1
        else:
            bedlist[1] = position - args.sgRNA_length - 1
            bedlist[2] = position - 1

    bedlist[1] = str(bedlist[1])
    bedlist[2] = str(bedlist[2])

    return "\t".join(bedlist)


def write_results(final_targets: List[str], gRNA_output_path: str, args: argparse.Namespace) -> None:
    """
    Writes the processed sgRNA target results to a specified output file.
    """
    with open(gRNA_output_path, "a") as k:
        for line in final_targets:
            if not args.coords_as_active_site:
                k.write(reverse_cut_site_offset(line, args))
            else:
                k.write(line)


def check_files(locations_to_keep: Optional[List[str]], locations_to_discard: Optional[List[str]]) -> bool:
    if not locations_to_keep:
        locations_to_keep = []
    if not locations_to_discard:
        locations_to_discard = []
    gtf_extensions = ["gtf", "gff", "gff2", "gff3"]
    other_exensions = ["bed"]
    found_gtf = False
    for file in locations_to_keep + locations_to_discard:
        ext = file.split(".")[-1].lower()
        if ext not in gtf_extensions + other_exensions:
            raise ValueError(f"\n\t{file} of unknown type. Please pass GTF/GFF/BED.")

        if ext in gtf_extensions:
            found_gtf = True

    return found_gtf


def main(args: Optional[argparse.Namespace] = None) -> None:
    if args is None:
        args = parse_arguments()

    args.fasta, was_gzipped = decompress_gzip_if_needed(args.fasta)

    # if args.active_site_offset_5 >= args.active_site_offset_3:
    #     raise ValueError("--active_site_offset_5 should be less than or equal to --active_site_offset_3")
    args.gc_range = sorted(args.gc_range)

    gRNA_output_path, _ = create_output(args.fasta, outdir=args.output_directory, extension="gRNA")
    gRNA_output_path += ".bed"

    with open(gRNA_output_path, "w") as f:
        f.write("#chr\tstart\tstop\tid,sequence,pam,chromosome,position,sense\tcontext\tstrand\n")

    keep_chroms = None
    if not args.locations_to_keep:
        print(
            "\n\tWARNING: '--locations_to_keep' flag not set.\n \
              All sgRNA across all fasta entries will be returned (except --locations_to_discard).\n \
              This can lead to a very large intermediate file being saved to disk.\n"
        )

    gtf = check_files(args.locations_to_keep, args.locations_to_discard)
    if args.feature and not gtf:
        print(f"\n\tNo GTF file, '--feature {args.feature}' will be ignored.\n")

    targets_to_keep = merge_targets(args.locations_to_keep, gtf_feature=args.feature, operation=args.join_operation)
    targets_to_discard = merge_targets(args.locations_to_discard, gtf_feature=args.feature, operation="merge")

    if targets_to_keep and targets_to_discard:
        targets_to_keep = targets_to_keep.subtract(targets_to_discard)
        targets_to_discard = None

    # if targets_to_keep:
    #     chroms = unique_chromosomes(targets_to_keep)

    if targets_to_keep:
        keep_chroms = get_chromosome_boundaries(targets_to_keep)
        print("\n\tChromosomes for which to find targets:\t" + " ".join(keep_chroms.keys()))

    pam_set = map_ambiguous_sequence(args.pam)
    rev_pam_set = list(map(revcom, pam_set))

    pad = args.sgRNA_length + args.context_window[0] + args.context_window[1] + max(len(p) for p in pam_set) + 10
    has_locations = bool(targets_to_keep or targets_to_discard)
    # With no region filtering each worker applies the cut-site -> protospacer transform
    # itself (in parallel). With locations, rows must stay in cut-site coords until after
    # intersect/subtract, then write_results applies the transform (matching the original).
    apply_reverse = (not has_locations) and (not args.coords_as_active_site)

    # Flat task list: tile every (kept) chromosome into windows, across ALL chromosomes.
    fa = pysam.FastaFile(args.fasta)
    tasks = []
    for chrom in fa.references:
        chrom_len = fa.get_reference_length(chrom)
        if keep_chroms and chrom not in keep_chroms:
            continue
        if chrom_len < args.min_chr_length:
            continue
        if keep_chroms:
            r_start, r_end = keep_chroms[chrom]
            r_start = max(r_start - max(args.context_window), 0)
            r_end = min(r_end + max(args.context_window), chrom_len)
        else:
            r_start, r_end = 0, chrom_len
        for w_start in range(r_start, r_end, args.chunk_size):
            tasks.append((chrom, w_start, min(w_start + args.chunk_size, r_end)))
    fa.close()
    print(f"\n\tScanning {len(tasks)} windows across {args.threads} workers\n")

    out_dir = os.path.dirname(gRNA_output_path) or "."
    tmpdir = tempfile.mkdtemp(prefix="ggwin_", dir=out_dir)
    paths = []
    with ProcessPoolExecutor(max_workers=args.threads, initializer=_init_worker, initargs=(args.fasta,)) as executor:
        futures = [
            executor.submit(process_window, c, ws, we, args, pam_set, rev_pam_set, pad, tmpdir, apply_reverse)
            for (c, ws, we) in tasks
        ]
        for future in futures:
            p = future.result()
            if p:
                paths.append(p)

    if not has_locations:
        # workers already wrote final rows (protospacer coords unless --coords_as_active_site)
        with open(gRNA_output_path, "ab") as out:
            for p in paths:
                with open(p, "rb") as ph:
                    shutil.copyfileobj(ph, out)
                os.remove(p)
    else:
        combined = os.path.join(tmpdir, "_combined.bed")
        with open(combined, "wb") as cf:
            for p in paths:
                with open(p, "rb") as ph:
                    shutil.copyfileobj(ph, cf)
                os.remove(p)
        bed = BedTool(combined)
        if targets_to_keep:
            final_targets = bed.intersect(targets_to_keep, u=True)
        else:
            final_targets = bed.subtract(targets_to_discard, A=True)
        write_results([str(line) for line in final_targets], gRNA_output_path, args)
        os.remove(combined)

    try:
        os.rmdir(tmpdir)
    except OSError:
        pass

    print(f"\n\tSaved output file to {gRNA_output_path}\n")

    if was_gzipped:
        remove_file(args.fasta)


if __name__ == "__main__":
    main()
