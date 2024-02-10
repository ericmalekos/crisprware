import pandas as pd
from pybedtools import BedTool


def adjust_interval_coordinates(interval, subtract_amount, add_amount, chrom_length=float('inf')):
    """
    Adjusts the start and end coordinates of an interval.

    Parameters:
    - interval (Bedtool interval)
    - subtract_amount (int): The amount to subtract from the interval's start coordinate.
    - add_amount (int): The amount to add to the interval's end coordinate.
    - chrom_length (int): length of chromosome to keep interval in bounds

    Returns:
    - Bedtool interval: The adjusted interval. Note that the original interval object is modified 
      and returned.
    """
    if subtract_amount < 0 or add_amount < 0:
        raise ValueError("subtract_amount and add_amount must be non-negative")
    
    interval.start = max(0, interval.start - subtract_amount)
    interval.end = min(interval.end + add_amount, chrom_length)
    return interval


def preprocess_file(file, gtf_feature="exon"):
    """
    Preprocesses a file and converts it to a sorted BedTool object.

    Parameters:
    - file (str): The path to the input GTF, GFF, or BED file.
    - gtf_feature (str): The feature type to filter GTF/GFF files by. 
      Defaults to 'exon'.

    Returns:
    - BedTool: A BedTool object representing the sorted genomic intervals from the input file.
    """
    extensions = ['.gtf', '.gff', '.gff2', '.gff3']
    extension = file[-4:].lower()
    if extension in extensions:
        df = pd.read_csv(file, sep='\t', header=None, comment='#', quoting=3)
        df = df[df[2] == gtf_feature]
        
        # Check if the DataFrame is empty after filtering
        if df.empty:
            raise ValueError(f"No features of type '{gtf_feature}' found in the file.")
        
        return BedTool.from_dataframe(df, na_rep=".", quoting=3).sort()
    else:
        return BedTool(file).sort()

def merge_targets(files, gtf_feature="exon", operation="intersect", window = [0,0]):
    """
    Intersect or merge an arbitrary number of GTF/GFF/BED files using pybedtools.

    Parameters:
    - files (list): List of paths to GTF/GFF/BED files to merge.
    - gtf_feature (str): The feature type to filter on in GTF/GFF files.
    - operation (str): The operation to perform ('intersect' or 'merge').
    - window (list): values to expand the interval by

    Returns:
    - BedTool: A BedTool object representing the merge of all files.
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
        
    result = result.merge()
    result = result.sort()

    return result

def filter_gtf_by_transcript_ids(input_file, transcript_ids):
    """
    Filters the GTF/GFF file based on the provided set of transcript IDs.

    Parameters:
    - input_file: path to the GTF/GFF file
    - transcript_ids: a set containing transcript IDs to be filtered for.

    Returns:
    - A list containing lines from the GTF/GFF file that match the given transcript IDs.
    """
    filtered_lines = []

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue

            _, attributes = parse_line(line)
            transcript_id = attributes.get('transcript_id', None)


            if transcript_id and transcript_id in transcript_ids:
                filtered_lines.append(line)

    return filtered_lines

def gtf_to_tss_bed(input_gtf, upstream=500, downstream=500):
    """
    Converts GTF file to a BED file with TSS positions for each transcript.

    Parameters:
    - input_gtf: Path to the input GTF file.
    - upstream: The number of bases upstream of the TSS to include in the BED entry.
    - downstream: The number of bases downstream of the TSS to include in the BED entry.

    Returns:
    - A list strings, each string is a line of a bed file
    """

    tss_entries = []

    with open(input_gtf, 'r') as gtf:
        for line in gtf:
            if not line.startswith("#"):  # skip header lines
                fields, attributes = parse_line(line)

                if fields[2] == 'transcript':
                    chrom = fields[0]
                    start = int(fields[3]) # 0-based start for BED format
                    end = int(fields[4])  # 1-based end for BED format
                    strand = fields[6]
                    
                    # Ensure the transcript_id attribute is present before proceeding
                    if 'transcript_id' not in attributes:
                        continue
                    
                    tss_start, window_start,window_end = -1,-1,-1
                    # Determine TSS and create a window around it based on strand orientation
                    if strand == '+':
                        window_start = max(0, start - upstream - 1)
                        window_end = start + downstream
                        tss_start = start - 1
                    else:
                        window_start = max(0, end - downstream - 1)
                        window_end = end + upstream
                        tss_start = end - 1

                    # Create BED entry
                    bed_entry = [chrom, str(window_start), str(window_end), attributes['transcript_id'], '0', strand,
                                 str(tss_start), str(tss_start+1)]
                    tss_entries.append('\t'.join(bed_entry))
    
    if not tss_entries:
        print('\tTSS could not be inferred.\n\tCheck that ' + input_gtf + ' has \'transcript\' attributes in the third column.')

    return tss_entries



def parse_input(input_file):
    """
    Parses a genomic data file and extracts gene-related information.

    Parameters:
    - input_file (str): The path to the input file containing genomic data. 
                        The file is expected to have specific fields in each line, 
                        representing genomic features like exons and CDS.

    Returns:
    - dict: A dictionary with gene_ids as keys. Each key maps to a dictionary that contains:
        - 'exons': A dictionary with transcript_ids as keys and a list of exon start-end tuples as values.
        - 'CDS_coords': A dictionary similar to 'exons', but for CDS coordinates.
        - 'CDS': A set of transcript_ids that have CDS.
        - 'strand': The strand ('+' or '-') of the gene.
        - 'chromosome': The chromosome on which the gene is located.
    """

    genes = {}

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields, attributes = parse_line(line)
            if not fields or not attributes:
                continue

            gene_id = attributes.get('gene_id', None)
            transcript_id = attributes.get('transcript_id', None)
            chromosome = fields[0]
            start, end = int(fields[3]), int(fields[4])
            feature_type = fields[2]
            strand = fields[6]

            if gene_id not in genes:
                genes[gene_id] = {'exons': {}, 'CDS_coords': {}, 'CDS': set(), 'strand': strand, 'chromosome': chromosome}

            if feature_type == 'exon':
                if transcript_id not in genes[gene_id]['exons']:
                    genes[gene_id]['exons'][transcript_id] = []
                genes[gene_id]['exons'][transcript_id].append((start, end))
            elif feature_type == 'CDS':
                genes[gene_id]['CDS'].add(transcript_id)
                if transcript_id not in genes[gene_id]['CDS_coords']:
                    genes[gene_id]['CDS_coords'][transcript_id] = []
                genes[gene_id]['CDS_coords'][transcript_id].append((start, end))

    return genes

def parse_line(line):
    """
    Parses a single line from a GTF/GFF file and extracts its fields and attributes.
    The function differentiates between GTF and GFF formats by checking for space (' ') or equal sign ('=') 
    delimiters in the attributes part of the line.

    Parameters:
    - line (str): A single line from a GTF/GFF file.

    Returns:
    - tuple: A tuple containing two elements:
        1. fields (list): A list of fields extracted from the line, typically including 
           chromosome, source, feature type, start, end, score, strand, and frame.
        2. attributes (dict): A dictionary of attributes parsed from the line. 
           Keys are attribute names (like 'gene_id' or 'transcript_id'), and values are the corresponding attribute values.
    """

    fields = line.strip().split('\t')
    attributes = {}

    for attr in fields[8].split(';'):
        attr = attr.strip()
        if not attr:
            continue

        # Checking if it's more like a GTF (space delimiter) or GFF (equal delimiter)
        if ' ' in attr:
            key, value = attr.split(' ', 1)
            attributes[key] = value.strip('"')
        elif '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key] = value

    return fields, attributes

def extract_transcript_gene_relationship(input_file):
    """
    Extracts a mapping between transcript_ids and gene_ids from the GTF/GFF file.

    Parameters:
    - input_file: path to the GTF/GFF file

    Returns:
    - A dictionary with transcript_ids as keys and gene_ids as values.
    """
    relationship = {}
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue

            _, attributes = parse_line(line)
            gene_id = attributes.get('gene_id', None)
            transcript_id = attributes.get('transcript_id', None)
            if transcript_id and gene_id:
                relationship[transcript_id] = gene_id
    return relationship

def write_utr(exon, start, end, cds_start, cds_end, strand, attributes):
    """
    Generates UTR (Untranslated Region) entries for a given exon based on the CDS (Coding Sequence) coordinates.

    Parameters:
    - exon (list): A list of fields representing an exon. Typically includes chromosome, source, feature type, 
      start, end, score, strand, frame, and attributes.
    - start (int): The start position of the exon.
    - end (int): The end position of the exon.
    - cds_start (int): The start position of the CDS.
    - cds_end (int): The end position of the CDS.
    - strand (str): The strand of the exon ('+' or '-').
    - attributes (str): A string of attributes for the exon.

    Returns:
    - str: A string representation of the UTR regions derived from the exon. Each UTR entry is in the same format as 
      the input exon fields, separated by tabs, and each entry is on a new line. The function can return entries for 
      5'UTR, 3'UTR, both, or none, depending on the exon and CDS positions.

    The function creates UTR entries by adjusting the start and end positions of the exon based on the overlap 
    with the CDS region. For the positive strand, it considers the region before the CDS start as 5'UTR and the region 
    after the CDS end as 3'UTR. For the negative strand, it considers the region after the CDS end as 5'UTR and the region 
    before the CDS start as 3'UTR.
    """
    utr_str = ""
    if strand == '+':
        if start < cds_start:
            utr5_fields = exon.copy()
            utr5_fields[2] = '5UTR'
            utr5_fields[4] = str(min(end, cds_start - 1))
            utr5_fields[8] = attributes
            utr_str += "\t".join(utr5_fields) + "\n"
        
        if end > cds_end:
            utr3_fields = exon.copy()
            utr3_fields[2] = '3UTR'
            utr3_fields[3] = str(max(start, cds_end + 1))
            utr3_fields[8] = attributes
            utr_str += "\t".join(utr3_fields) + "\n"
    else:
        if end > cds_end:
            utr5_fields = exon.copy()
            utr5_fields[2] = '5UTR'
            utr5_fields[3] = str(max(start, cds_end + 1))
            utr5_fields[8] = attributes
            utr_str += "\t".join(utr5_fields) + "\n"
        
        if start < cds_start:
            utr3_fields = exon.copy()
            utr3_fields[2] = '3UTR'
            utr3_fields[4] = str(min(end, cds_start - 1))
            utr3_fields[8] = attributes
            utr_str += "\t".join(utr3_fields) + "\n"

    return utr_str

def generate_output_str(genes, consensus_exons, consensus_CDS, label = "consensus" ):
    """
    Generates a string representation of consensus gene models in GTF format.
    Called by create_constitutive_model and create_metagene_model

    Parameters:
    - genes (dict): A dictionary containing gene information, including chromosome, strand, and potentially other metadata.
                    The keys are gene IDs, and values are dictionaries with details such as chromosome and strand.
    - consensus_exons (dict): A dictionary with gene IDs as keys and sets of exon start-end tuples as values, 
                              representing consensus exon regions for each gene.
    - consensus_CDS (dict): A dictionary similar to consensus_exons but for coding sequences. If a gene has consensus CDS,
                            it will be used to define CDS regions and calculate UTRs.
    - label (str, optional): A label to annotate the consensus gene models. Default is "consensus".

    Returns:
    - str: A string containing the formatted GTF lines for the consensus gene models. This includes one line for each 
           transcript, exon, and CDS, with UTRs included if consensus CDS data is available.
    """

    output_str = ''

    for gene_id, exons in consensus_exons.items():
        # Define CDS bounds based on consensus_CDS
        cds_list = sorted(list(consensus_CDS.get(gene_id, [])))
        if cds_list:
            cds_start, _ = cds_list[0]
            _, cds_end = cds_list[-1]
        else:
            cds_start, cds_end = None, None

        if exons:
            exons = sorted(list(exons), key=lambda x: x[0])
            transcript_start = exons[0][0]
            transcript_end = exons[-1][1]
            transcript_line = [
                genes[gene_id]['chromosome'], label, 'transcript', str(transcript_start), str(transcript_end), '.', 
                genes[gene_id]['strand'], '.', f'gene_id "{gene_id}"; transcript_id "{gene_id}.{label}"'
            ]
            output_str += "\t".join(transcript_line) + "\n"

        for exon_tuple in exons:
            start, end = exon_tuple
            exon_line = [
                genes[gene_id]['chromosome'], label, 'exon', str(start), str(end), '.', 
                genes[gene_id]['strand'], '.', f'gene_id "{gene_id}"; transcript_id "{gene_id}.{label}"'
            ]
            output_str += "\t".join(exon_line) + "\n"
            if cds_start and cds_end:  # Only write UTRs if there is a consensus CDS
                attributes = f'gene_id "{gene_id}"; transcript_id "{gene_id}.{label}"'
                output_str += write_utr(exon_line, start, end, cds_start, cds_end, genes[gene_id]['strand'], attributes)

        for cds_tuple in consensus_CDS.get(gene_id, []):
            start, end = cds_tuple
            cds_line = [
                genes[gene_id]['chromosome'], label, 'CDS', str(start), str(end), '.', 
                genes[gene_id]['strand'], '.', f'gene_id "{gene_id}"; transcript_id "{gene_id}.{label}"'
            ]
            output_str += "\t".join(cds_line) + "\n"

    return output_str.strip()



def is_inside(region1, region2):
    """
    Check if region1 is inside region2.
    """
    return region2[0] <= region1[0] and region2[1] >= region1[1]




#TODO CLEANUP create_constitutive_model & create_metagene_model functions. Messy and redundant currently

def create_constitutive_model(input_file):
    """
    Generates a constitutive model of gene expression from a GTF file.

    This function processes a GTF file to identify consensus exons and coding sequences (CDS) across all transcripts of a gene.
    It aims to construct a model where only those exons and CDS that are present in all transcripts (or in all transcripts containing a CDS)
    are considered as constitutive, i.e., consistently expressed across different isoforms of the gene.

    The process involves two main phases:
    1. Parsing the input GTF file to extract gene, exon, and CDS information into a structured format.
    2. Iterating over this structured data to identify consensus regions across transcripts for each gene.

    Parameters:
    - input_file (str): Path to the input GTF file containing annotations for genes, transcripts, exons, and CDS.

    Returns:
    - tuple: A tuple containing two elements:
        1. output_str (str): A string representation of the consensus gene model, formatted similarly to GTF, including
           consensus exons and CDS for each gene.
        2. genes_without_consensus (set): A set of gene IDs for which no consensus could be established, indicating
           variability in exon or CDS presence across transcripts.
    """
    genes = parse_input(input_file)

    # Phase 2: Find consensus exons and CDS
    consensus_exons = {}
    consensus_CDS = {}

    for gene_id, data in genes.items():
        #cds_transcripts = get_transcripts_with_cds(data['CDS'])
        cds_transcripts = set(data['CDS'])
        all_transcripts = set(data['exons'].keys())
        
        if not cds_transcripts:  # If no transcripts have CDS, continue to next gene
            for transcript, exons in data['exons'].items():
                for exon in exons:
                    consensus = overlapping_regions_for_transcripts(exon, all_transcripts, data['exons'])
                    #print(consensus)
                    if consensus:
                        consensus_exons[gene_id] = consensus_exons.get(gene_id, set())
                        consensus_exons[gene_id].add(consensus)
            
        # Process exons
        else:
            for transcript, exons in data['exons'].items():
                for exon in exons:
                    consensus = overlapping_regions_for_transcripts(exon, cds_transcripts, data['exons'])
                    #print(consensus)
                    if consensus:
                        consensus_exons[gene_id] = consensus_exons.get(gene_id, set())
                        consensus_exons[gene_id].add(consensus)

            # Process CDS
            for transcript, cds_coords in data['CDS_coords'].items():
                for cds in cds_coords:
                    cds_consensus = overlapping_regions_for_transcripts(cds, cds_transcripts, data['CDS_coords'])
                    if cds_consensus:
                        consensus_CDS[gene_id] = consensus_CDS.get(gene_id, set())
                        consensus_CDS[gene_id].add(cds_consensus)

    output_str = generate_output_str(genes, consensus_exons, consensus_CDS)
    genes_without_consensus = set(genes.keys()) - set(consensus_exons.keys())

    return output_str, genes_without_consensus
        
def create_metagene_model(input_file):
    """
    Generates a metagene model from gene annotations in a GTF file.

    This function processes a GTF file to create a metagene model for each gene. The metagene model combines all
    exon and coding sequence (CDS) regions across all transcripts of a gene into merged regions. This approach
    simplifies the gene structure by creating a unified representation of exons and CDS that appear across
    different transcripts, facilitating analyses that require a generalized view of gene architecture.

    The metagene model construction involves the following steps:
    1. Parsing the input GTF file to extract detailed information about genes, including their exons and CDS.
    2. For each gene, combining exon regions from all (or CDS-containing) transcripts and merging overlapping
       or adjacent regions to form a unified set of exons for the metagene model.
    3. Similarly, if applicable, combining and merging CDS regions from transcripts that contain CDS.
    4. Generating a string representation of the metagene model in GTF format.

    Parameters:
    - input_file (str): Path to the GTF file containing the gene annotations to be processed.

    Returns:
    - str: A string containing the metagene model in GTF format. This includes entries for the merged
           exons and CDS of each gene, annotated with a 'metagene' label.
    """
    
    genes = parse_input(input_file)
    
    # Phase 2: Combine all exons for metagene model
    metagene_exons = {}
    metagene_CDS = {}

    for gene_id, data in genes.items():
 
        #cds_transcripts = get_transcripts_with_cds(data['CDS'])
        cds_transcripts = set(data['CDS'])

        all_transcripts = set(data['exons'].keys())

        if len(all_transcripts) == 0:
            print("\t\tWARNING:\tNo exon found for " + gene_id + ", skipping this gene. Check GTF.")
            continue

        use_transcripts = cds_transcripts if cds_transcripts else all_transcripts

        combined_exons = []

        for transcript in use_transcripts:
            combined_exons.extend(data['exons'][transcript])
              
        combined_exons.sort(key=lambda x: x[0])
        
        merged_exons = []
        current_exon = combined_exons[0]

        for exon in combined_exons[1:]:
            if exon[0] <= current_exon[1] + 1:
                current_exon = (current_exon[0], max(current_exon[1], exon[1]))
            else:
                merged_exons.append(current_exon)
                current_exon = exon
        merged_exons.append(current_exon)

        metagene_exons[gene_id] = merged_exons

        # If CDS transcripts are present, combine and merge them similarly to exons
        if cds_transcripts:
            combined_cds = []
            for transcript in cds_transcripts:
                combined_cds.extend(data['CDS_coords'][transcript])
            combined_cds.sort(key=lambda x: x[0])

            merged_cds = []
            current_cds = combined_cds[0]
            for cds in combined_cds[1:]:
                if cds[0] <= current_cds[1] + 1:
                    current_cds = (current_cds[0], max(current_cds[1], cds[1]))
                else:
                    merged_cds.append(current_cds)
                    current_cds = cds
            merged_cds.append(current_cds)
            metagene_CDS[gene_id] = merged_cds

    # Phase 3: save output (remains mostly unchanged, but work with `metagene_exons` and `metagene_CDS` dicts)
    output_str = generate_output_str(genes, metagene_exons, metagene_CDS, label='metagene')

    return output_str


def overlapping_regions_for_transcripts(region, transcripts, all_regions):
    """
    Identifies the common overlapping region across specified transcripts.
    Helper function for create_constitutive_model.

    Parameters:
    - region (tuple): A tuple specifying the start and end positions of the reference region (start, end).
    - transcripts (iterable): An iterable of transcript identifiers to be considered for finding the overlap.
    - all_regions (dict): A dictionary mapping transcript identifiers to lists of regions (tuples of start and end positions)
                          associated with that transcript.

    Returns:
    - tuple or None: If a common overlapping region exists, returns a tuple (start, end) representing the maximal 
                     start position and minimal end position common to all specified transcripts. If no common 
                     overlapping region is found, returns None.

    """

    overlapping_regions = [region]
    for transcript in transcripts:
        overlaps_current_transcript = False
        for t_region in all_regions.get(transcript, []):
            if t_region[1] >= region[0] and t_region[0] <= region[1]:
                overlapping_regions.append(t_region)
                overlaps_current_transcript = True
        if not overlaps_current_transcript:
            return None
    return get_max_start_min_end(overlapping_regions)


def get_max_start_min_end(regions):
    """
    Determines the common overlapping region across a list of regions.

    Parameters:
    - regions (list of tuples): A list where each tuple contains two integers representing the start and 
                                end positions of a region (start, end).

    Returns:
    - tuple or None: If a common overlapping region is found, returns a tuple (max_start, min_end) representing 
                     this region. If no common overlap exists (indicating disjoint regions), returns None.
    """
    if not regions:
        return None
    max_start = max([start for start, _ in regions])
    min_end = min([end for _, end in regions])
    
    if max_start > min_end:
        return None
    
    return (max_start, min_end)

