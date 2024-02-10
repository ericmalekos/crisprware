from pathlib import Path
import pandas as pd
from pybedtools import BedTool
from itertools import product
from Bio import SeqIO

# Constants
NTS = list("ACTG")
NT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}
AMBIGUITY_CODES = {
    "R": "AG", "Y": "CT", "K": "GT", "M": "AC", 
    "S": "GC", "W": "AT", "B": "CGT", "D": "AGT", 
    "H": "ACT", "V": "ACG", "N": "ACGT"
}


def map_ambiguous_sequence(sequence):
    """
    Map an ambiguous DNA sequence to all possible non-ambiguous sequences.

    Parameters:
    sequence (str): The ambiguous DNA sequence.

    Returns:
    list of str: A list of all possible non-ambiguous sequences.
    """

    # Generate a list of lists where each ambiguous ntleotide is replaced by a list of possible nucleotides
    possible_nucleotides = [[nt] if nt not in AMBIGUITY_CODES else list(AMBIGUITY_CODES[nt]) for nt in sequence.upper()]

    # Get all possible sequences by calculating the cartesian product of the lists of possible nucleotides
    possible_sequences = [''.join(seq) for seq in product(*possible_nucleotides)]

    return possible_sequences


def revcom(dna):
    return "".join(list(map(lambda n: NT_MAP[n], list(dna)))[::-1])

def preprocess_file(file, gtf_feature="exon"):
    extensions = ['.gtf', '.gff', '.gff2', '.gff3']
    extension = Path(file).suffix.lower()
    if extension in extensions:
        df = pd.read_csv(file, sep='\t', header=None, comment='#', quoting=3)
        df = df[df[2] == gtf_feature]
        return BedTool.from_dataframe(df,na_rep=".",quoting=3)
    else:
        return BedTool(file)

def merge_targets(files, gtf_feature="exon", operation="intersect", window = [0,0]):
    """
    Intersect or merge an arbitrary number of GTF/GFF/BED files using pybedtools.

    Parameters:
    files (list): List of paths to GTF/GFF/BED files to merge.
    gtf_feature (str): The feature type to filter on in GTF/GFF files.
    operation (str): The operation to perform ('intersect' or 'merge').

    Returns:
    BedTool: A BedTool object representing the merge of all files.
    """

    if not files:
        return None

    # Initialize the result with the first preprocessed BED file
    result = preprocess_file(files[0], gtf_feature)

    # Iteratively intersect the result with each subsequent preprocessed BED file
    for file in files[1:]:
        if operation == "intersect":
            result = result.intersect(preprocess_file(file=file, gtf_feature=gtf_feature))
        elif operation == "merge":
            result = result.cat(preprocess_file(file=file, gtf_feature=gtf_feature), postmerge=True)
        else:
            raise ValueError("Invalid operation. operation should be either 'merge' or 'intersect'.")

    if window[0] != 0 or window[1] != 0:
        result = result.each(adjust_interval_coordinates, 
                    subtract_amount=window[0], add_amount=window[1])

    result = result.sort()
    result = result.merge()

    return result


def subset_fasta_with_bed(fasta_path, bed_path, output_fasta_path):
    """
    Subsets a FASTA file using regions specified in a BED file.
    Handles cases where the BED file has chromosomes not in the FASTA,
    or if the end of the BED interval is beyond the length of the FASTA sequence.
    The output FASTA record names include the coordinates used in subsetting.

    Parameters:
    fasta_path (str): Path to the input FASTA file.
    bed_path (str): Path to the BED file specifying regions.
    output_fasta_path (str): Path to output the subsetted FASTA file.
    """

    # Read BED file using pybedtools
    bed = BedTool(bed_path)

    # Read the FASTA file
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'))

    # Create a new FASTA file for output

    with open(output_fasta_path, 'w') as output_fasta:
        for interval in bed:
            seq_id = interval.chrom
            start = int(interval.start)
            end = int(interval.end)
            
            # Check if the sequence ID is in the FASTA file
            if seq_id in fasta_sequences:
                seq_length = len(fasta_sequences[seq_id].seq)

                # Adjust the end position if it is beyond the sequence length
                if end > seq_length:
                    print(f"\n\tEnd coordinate for interval {seq_id} {start} {end}\n\tis beyond sequence length of {seq_length}\n\tReducing end coordinate")
                    end = seq_length

                # Extract the sequence from the FASTA
                sequence = fasta_sequences[seq_id].seq[start:end]

                # Format the record name to include coordinates
                record_id = f"{seq_id}_{start}:{end}"

                SeqIO.write(SeqIO.SeqRecord(sequence, id=record_id, description=""), output_fasta, "fasta-2line")
            else:
                print(f"\n\tWarning: Chromosome {seq_id} not found in the FASTA file.\n")

def remove_restricted(sgRNA, patterns, flank_5, flank_3):
    """Calculate the GC content of a nucleotide sequence."""
    context = flank_5 + sgRNA + flank_3
    for pattern in patterns:
        if (pattern in context) or (pattern in revcom(context)):
            return True

    return False

def calculate_gc_content(sequence):
    """Calculate the GC content of a nucleotide sequence."""
    
    gc_count = sequence.count('G') + sequence.count('C')
    
    total_count = len(sequence)
    
    # Calculate the GC content
    gc_content = (gc_count / total_count) * 100
    
    return gc_content

def include_sgRNA(args, sgRNA):
    """
    TODO: pass arguments individually instead of as "args"
    
    Determines whether a given sgRNA sequence should be included based on GC content,
    presence of poly-T sequences, and restriction patterns.

    Parameters:
    - args (Namespace): A Namespace object containing the following attributes:
        - gc_range (tuple): A tuple specifying the minimum and maximum allowable GC content
          percentages for the sgRNA sequence.
        - discard_poly_T (bool): If True, sgRNA sequences containing the sequence 'TTTT'
          will be excluded.
        - restriction_patterns (list): A list of strings representing the recognition
          sequences of restriction enzymes. sgRNA sequences containing any of these
          patterns within specified flanking regions will be excluded.
        - flank_5 (int): The length of the 5' flanking region to consider for restriction
          pattern presence.
        - flank_3 (int): The length of the 3' flanking region to consider for restriction
          pattern presence.
    - sgRNA (dict): A dictionary representing a single sgRNA, with at least the key
        'sequence' containing the sgRNA nucleotide sequence as a string.

    Returns:
    - bool: True if the sgRNA meets all inclusion criteria, False otherwise.
    """

    gc = calculate_gc_content(sgRNA["sequence"])

    if gc < args.gc_range[0] or gc > args.gc_range[1]:
        return False

    if args.discard_poly_T and "TTTT" in sgRNA["sequence"]:
        return False

    if len(args.restriction_patterns) > 0 and \
        remove_restricted(sgRNA["sequence"], args.restriction_patterns, \
                            args.flank_5, args.flank_3):
        return False    
    
    return True

def get_chromosome_boundaries(bed_obj):
    """x
    Returns a dictionary with chromosomes as keys and a tuple (smallest_start, largest_end) as values.

    Parameters:
    - bed_obj (pybedtools.BedTool): A BedTool object containing intervals.

    Returns:
    - dict: Chromosomes as keys and boundaries as values.
    """
    
    # Sort the BedTool object
    sorted_bed = bed_obj.sort()

    # Dictionary to store results
    boundaries = {}

    # Variables to store current chromosome data
    current_chrom = None
    smallest_start = None
    largest_end = None

    # Iterate over intervals in the sorted BedTool object
    for interval in sorted_bed:
        if current_chrom is None:
            # First chromosome being processed
            current_chrom = interval.chrom
            smallest_start = interval.start
            largest_end = interval.end
        elif current_chrom == interval.chrom:
            # Update the largest end for the current chromosome
            largest_end = max(largest_end, interval.end)
        else:
            # Save results for the completed chromosome
            boundaries[current_chrom] = (smallest_start, largest_end)

            # Start processing the next chromosome
            current_chrom = interval.chrom
            smallest_start = interval.start
            largest_end = interval.end

    # Save results for the last chromosome
    boundaries[current_chrom] = (smallest_start, largest_end)

    return boundaries

def unique_chromosomes(input):
    """
    Return a set of unique chromosome IDs from a BED, GTF, or GFF file.
    """

    chroms = set()

    if isinstance(input, BedTool):
        chroms = {interval.chrom for interval in input}
    else:
        with open(input, 'r') as file:
            for line in file:
                # Skip comments, especially relevant for GTF and GFF
                if line.startswith("#"):
                    continue

                fields = line.strip().split('\t')
                chroms.add(fields[0])

    return chroms