#!/usr/bin/env python3

'''
    This modifies "generate_kmers.py" script from the Guidescan repo

    EX Usage:
        generate_guides --pam NGG  --sgRNA_length 20  -f ./chr21Index/chr21.fa --discard_poly_T --threads 8

        generate_guides --sgRNA_length 20  -f ./chr21Index/chr21.fa --discard_poly_T --threads 8 --locations_to_keep ./AnnotationFiles/NCAM2.gtf -o "NCAM_"


'''
from Bio import SeqIO
import argparse
from concurrent.futures import ProcessPoolExecutor
from pybedtools import BedTool
from utils.dna_sequence_functions import NTS, map_ambiguous_sequence, \
    revcom, merge_targets, include_sgRNA, get_chromosome_boundaries
from utils.utility_functions import create_output, decompress_gzip_if_needed, remove_file

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generates a set of sgRNA matching user specifications. \
            Default settings are for active SpCas9 [--pam=NGG --active_site_offset_5=-4 --active_site_offset_3=-3 --sgRNA_length=20]"
    )

    parser.add_argument(
        "-f", "--fasta",
        type=str,
        help="FASTA file to use as a reference for sgRNA generation.",
        required=True
    )

    parser.add_argument(
        "-p", "--pam",
        type=str,
        help="Protospacer adjacent motif to match.\
             All IUPAC ambiguity codes are accepted as well as standard ATCG. [default: NGG]",
        default="NGG"
    )

    parser.add_argument(
        "-l", "--sgRNA_length",
        help="Length of sgRNA to generate. [default: 20]",
        type=int,
        default=20
    )


    parser.add_argument(
        "-w", "--context_window",
        nargs=2, 
        type=int, 
        help="Pass two, space-separated, integers to specifiy the nucleotide window around the sgRNA \
            as '<upstream>' '<downstream>'. This can be used for downstream scoring, \
            For Ruleset3 use -w 4 6 to obtain an appropriate score context. [default: 4 6]",
        default=[4,6]
    )

    parser.add_argument(
        "-5", "--active_site_offset_5",
        type=int,
        help="Where cut occurs relative to PAM 5' end. \
             To avoid error, use '=' sign when passing a negative number, \
                 e.g. --active_site_offset_5=-1 [default: -4]",
        default="-4"
    )

    parser.add_argument(
        "-3", "--active_site_offset_3",
        type=int,
        help="Where cut occurs relative to PAM 5' end. [default: -2] \
            To avoid error, use '=' sign when passing a negative number, \
            e.g. --active_site_offset_3=-3 [default: -4]",
        default="-4"
    )

    parser.add_argument(
        "-k", "--locations_to_keep",
        help="List of BED/GTF files with coordinates in \
            which the sgRNA desired. If the sgRNA cutsite does not intersect \
            coordinates in these files they are discarded. Leave blank \
            to keep all sgRNA. \
            e.g. atac_peak.bed genes.gtf",
        type=str,   
        default="",
        nargs='*'
    )

    parser.add_argument(
        "--feature",
        type=str,
        help="For any GTF/GFF in '--locations_to_keep', only this \
            feature will be used for determining appropriate sgRNA. \
            The feature should match an entry in the third column of the GTF/GFF. \
            [default: 'exon']",
        default="exon"
    )
    
    parser.add_argument(
        "--join_operation",
        type=str, 
        choices=["merge", "intersect"], 
        help="How to treat '--locations_to_keep' if multiple files are passed. \
        Either 'merge' or 'intersect' can be used and work as described in Bedtools. \
        If 'merge',  sgRNA will be kept if its cutsite intersects an entry in ANY of the files, if \
        'intersect' the cutsite  must intersect an entry in EACH file. [default: 'intersect']",
        default="intersect"
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
        nargs='*'
    )

    parser.add_argument(
        "--prefix",
        help="Prefix to use for sgRNA identifiers. [default: None]",
        type=str,
        default=""
    )

    parser.add_argument(
        "--gc_range",
        nargs=2, 
        type=int, 
        help="Pass two, space-separated, integers to specifiy the percentile \
         range of GC content e.g. '--gc_range 25 75'. [default: 0 100]",
        default=[0,100]
    )

    parser.add_argument(
        "--discard_poly_T",
        help="Whether to discard polyT (>TTT) sgRNA. \
            Recommend True for PolIII promoters [default: False]",
        default=False,
        action="store_true"
    )

    parser.add_argument(
        "--discard_poly_G",
        help="Whether to discard polyT (>GGGG) sgRNA. [default: False]",
        default=False,
        action="store_true"
    )
    
    parser.add_argument(
        "--restriction_patterns",
        help="Reject sgRNA with these restriction patterns. \
            Also checks 5'flank+sgRNA+3'flank, and reverse complement, if provided. \
            For multiple values, separate by space. \
            e.g. GCGGCCGC TCTAGA CACCTGC",
        type=str,   
        default="",
        nargs='*'
    )

    parser.add_argument(
        "--flank_5",
        help="include the 5' context of the lentivirus vector. \
            Used in conjunction with --restriction_patterns to \
            remove incompatible sgRNA",
        type=str,
        default=""
    )
    parser.add_argument(
        "--flank_3",
        help="include the 3' context of the lentivirus vector. \
            Used in conjunction with --restriction_patterns to \
            remove incompatible sgRNA",
        type=str,
        default=""
    )

    parser.add_argument(
        "" "--min_chr_length",
        help="Minimum chromosome length to consider for sgRNA generation. [default: 20]",
        type=int,
        default=20
    )

    parser.add_argument(
        "--pam_5_prime",
        help="If the PAM is positioned 5' to the protospacer \
        set this flag, e.g. for Cas12a sgRNAs [default: False]",
        default=False,
        action="store_true"
    )


    parser.add_argument(
        "--coords_as_active_site",
        help="Whether to output bed coordinates at the active site \
        rather than the coordinates of the entire protospacer. For purposes \
        of keeping or discarding sgRNAs, overlap with the active site \
        coordinates will be used regardless [default: True]",
        default=True,
        action="store_false"
    )

    parser.add_argument(
        "-t", "--threads",
        type=int,
        help="Number of threads. [default: 4]",
        default=4,
    )

    parser.add_argument(
        "-o", "--output_directory",
        help="Path to output. [default: current directory]",
        type=str,
        default=""
    )

    return parser.parse_args()

def find_sgRNA(args, pam, chrm, start, end, forward=True):
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
                sgRNA = chrm[index - k:index]
                context = chrm[index - k - upstr:index + downstr]
                position = index - k + 1
            else:
                sgRNA = chrm[index + len(pam):index + k + len(pam)]
                context = chrm[index - (downstr-len(pam)):index + k + len(pam) + upstr]
                position = index + 1
        else:
            if forward:
                sgRNA = chrm[index + len(pam):index + len(pam) + k]
                context = chrm[index - upstr + len(pam):index + len(pam) + k + downstr]
                position = index + 1
            else:
                sgRNA = chrm[index - k:index]
                context = chrm[index - k - downstr:index + upstr]
                position = index + 1

        index += 1

        if position < 0:
            continue

        yield sgRNA.upper(), position, context.upper()


def output_bed_line(args, chrm_name, sgRNA):
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
    guideScanEntry = ",".join([identifier, str(sgRNA['sequence']),
           sgRNA['pam'], chrm_name,
           str(sgRNA['position']), sgRNA['sense']])
    
    # Find the position of the Cleavage site 
    if not args.pam_5_prime:    
        if sgRNA['sense'] == "+":
            cleavePos1 = sgRNA['position'] + sgRNA["length"] + args.active_site_offset_5 - 1
            cleavePos2 = sgRNA['position'] + sgRNA["length"] + args.active_site_offset_3 - 1
        else:
            cleavePos1 = sgRNA['position'] - args.active_site_offset_3 + len(args.pam) - 1
            cleavePos2 = sgRNA['position'] - args.active_site_offset_5 + len(args.pam) - 1
    else:
        if sgRNA['sense'] == "+":
            cleavePos1 = sgRNA['position'] + len(args.pam) + args.active_site_offset_5 - 1
            cleavePos2 = sgRNA['position'] + len(args.pam) + args.active_site_offset_3 - 1
        else:
            cleavePos1 = sgRNA['position'] - args.active_site_offset_3 - 1
            cleavePos2 = sgRNA['position'] - args.active_site_offset_5 - 1
        
    bedLine = [chrm_name, str(cleavePos1), str(cleavePos2), guideScanEntry, sgRNA['context'], sgRNA['sense']]
    return "\t".join(bedLine)

def process_pam(args, pam, record, start, end, pam_set, rev_pam_set):
    """
    Processes PAM sequences within a specified genomic region and generates results
    for sgRNA candidates that meet the defined criteria.

    Parameters:
    - args (Namespace): A Namespace object containing parameters for sgRNA and context
        extraction, including sgRNA length, context window size, and PAM orientation.
    - pam (str): The PAM sequence to search for.
    - record (SeqRecord): A BioPython SeqRecord object representing the genomic region
        to be searched. Contains sequence (`seq`) and identifier (`id`).
    - start (int): The start position within the genomic region to begin the search.
    - end (int): The end position within the genomic region to end the search.
    - pam_set (set): A set of valid PAM sequences for forward strand analysis.
    - rev_pam_set (set): A set of valid PAM sequences for reverse strand analysis.

    Returns:
    - list: A list of formatted results for each valid sgRNA candidate found. Each result
        includes sgRNA sequence, position, PAM, sense, length, and context information,
        formatted as specified by the `output_bed_line` function.

    Raises:
    - ValueError: If the provided PAM sequence is not found in either `pam_set` or
        `rev_pam_set`.

    """
    
    results = []

    #if not args.pam_5_prime:
    context_length = args.sgRNA_length + args.context_window[0] + args.context_window[1]
    # if the pam is 5' of the protospacer, the context length needs to be adjusted by length of the pam
    #else:
    #    context_length = args.sgRNA_length + args.context_window[0] + args.context_window[1] + len(pam)
    chrm_seq = str(record.seq).upper()
    chrm_name = record.id

    #Process pams
    if pam in pam_set:
        for sgRNA, pos, context in find_sgRNA(args=args, pam=pam, chrm=chrm_seq, start=start, end=end, forward=True):
            if len(sgRNA) != args.sgRNA_length or len(context) != context_length: continue
            if not all(nt in NTS for nt in context): continue
            if include_sgRNA(args, {"sequence" : sgRNA, "position" : pos, "pam" : pam,
                    "sense": "+", "length" : args.sgRNA_length, "context": context}):
                results.append(output_bed_line(args, chrm_name,
                                                {"sequence" : sgRNA, "position" : pos, "pam" : args.pam,
                                                "sense": "+", "length" : args.sgRNA_length, "context": context}))
    elif pam in rev_pam_set:
        for sgRNA, pos, context in find_sgRNA(args=args, pam=pam, chrm=chrm_seq, start=start, end=end, forward=False):
            if len(sgRNA) != args.sgRNA_length or len(context) != context_length: continue
            if not all(nt in NTS for nt in context): continue
            sgRNA = revcom(sgRNA)
            context = revcom(context)
            if include_sgRNA(args, {"sequence" : sgRNA, "position" : pos, "pam" : pam,
                    "sense": "-", "length" : args.sgRNA_length, "context": context}):
                results.append(output_bed_line(args, chrm_name, 
                                                {"sequence" : sgRNA, "position" : pos, "pam" : args.pam,
                                                "sense": "-", "length" : args.sgRNA_length, "context": context}))
    else:
        raise ValueError(f"The PAM sequence '{pam}' is not found in the list of valid PAMs.")

    return results

def reverse_cut_site_offset(bedline, args):
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

def write_results(final_targets, gRNA_output_path, args):
    """
    Writes the processed sgRNA target results to a specified output file.
    """
    with open(gRNA_output_path, 'a') as k:
        for line in final_targets:
            if not args.coords_as_active_site:
                k.write(reverse_cut_site_offset(line, args))
            else:
                k.write(line)


def check_files(locations_to_keep, locations_to_discard):
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



def main():
    args = parse_arguments()


    args.fasta, was_gzipped = decompress_gzip_if_needed(args.fasta)

    # if args.active_site_offset_5 >= args.active_site_offset_3:
    #     raise ValueError("--active_site_offset_5 should be less than or equal to --active_site_offset_3")
    args.gc_range = sorted(args.gc_range)

    gRNA_output_path, _ = create_output(args.fasta, outdir=args.output_directory, extension="gRNA")
    gRNA_output_path += ".bed"

    with open(gRNA_output_path, 'w') as f:
        f.write('#chr\tstart\tstop\tid,sequence,pam,chromosome,position,sense\tcontext\tstrand\n')

    keep_chroms = None
    if not args.locations_to_keep:
        print("\n\tWARNING: '--locations_to_keep' flag not set.\n \
              All sgRNA across all fasta entries will be returned (except --locations_to_discard).\n \
              This can lead to a very large intermediate file being saved to disk.\n")
    

    gtf = check_files(args.locations_to_keep, args.locations_to_discard)
    if args.feature and not gtf:
        print(f"\n\tNo GTF file, '--feature {args.feature}' will be ignored.\n")  


    targets_to_keep = merge_targets(args.locations_to_keep, gtf_feature=args.feature , operation = args.join_operation)
    targets_to_discard = merge_targets(args.locations_to_discard, gtf_feature=args.feature , operation = "merge")

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
    

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        for record in SeqIO.parse(args.fasta, "fasta"):
            chr_results = []

            #skip chromosomes
            if keep_chroms and record.id not in keep_chroms:
                continue

            print("\tProcessing " + record.id)
            #skip short chromosomes
            if len(record.seq) < args.min_chr_length:
                print("\tSkipping record, too short:\t" + record.id)
                continue
            if keep_chroms:
                start, end = keep_chroms[record.id]
                start -= max(args.context_window) 
                start = max(start, 0)
                end += max(args.context_window)
                end = min(end, len(record.seq))
            else:
                start = 0
                end = len(record.seq)


            # Submitting each PAM for processing as separate task
            futures = []

            for pam in (pam_set + rev_pam_set):
                futures.append(executor.submit(process_pam, args=args, pam=pam, record=record, start=start, end=end, 
                                                pam_set=pam_set, rev_pam_set=rev_pam_set))

            # Waiting for all futures to complete and collecting the results
            for future in futures:
                chr_results.extend(future.result())

            bed_chr_results = BedTool("\n".join(chr_results), from_string=True)

            if targets_to_keep:
                final_targets = bed_chr_results.intersect(BedTool(targets_to_keep), u=True)
            elif targets_to_discard:
                final_targets = bed_chr_results.subtract(BedTool(targets_to_discard), A=True)
            else:
                final_targets = bed_chr_results

            final_targets_str = [str(line) for line in final_targets]
            write_results(final_targets_str, gRNA_output_path, args)

    if was_gzipped:
        remove_file(args.fasta)

if __name__ == "__main__":
    main()