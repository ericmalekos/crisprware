import unittest
from utils.gtf_bed_processing_functions import preprocess_file,merge_targets,\
    gtf_to_tss_tes_bed,parse_line, extract_transcript_gene_relationship, parse_input,\
        create_metagene_model,create_constitutive_model, truncate_gtf
import pandas as pd
from io import StringIO

class TestGTFBEDFunctions(unittest.TestCase):
    
    def test_preprocess_file(self):
        gtf = "./tests/test_data/test.gtf"
        result = preprocess_file(gtf, gtf_feature="exon")
        expected_result =   'test_chr\tTest\texon\t20\t35\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";\n'\
                            'test_chr\tTest\texon\t40\t50\t.\t+\t.\tgene_id "gene1"; transcript_id "transcript1";\n'\
                            'test_chr\tTest\texon\t120\t140\t.\t+\t.\tgene_id "gene2"; transcript_id "transcript2";\n'\
                            'test_chr\tTest\texon\t150\t170\t.\t+\t.\tgene_id "gene2"; transcript_id "transcript2";\n'\
                            'test_chr\tTest\texon\t200\t230\t.\t-\t.\tgene_id "gene3"; transcript_id "transcript3";\n'\
                            'test_chr\tTest\texon\t200\t230\t.\t-\t.\tgene_id "gene3"; transcript_id "transcript4";\n'\
                            'test_chr\tTest\texon\t240\t260\t.\t-\t.\tgene_id "gene3"; transcript_id "transcript3";\n'\
                            'test_chr\tTest\texon\t240\t260\t.\t-\t.\tgene_id "gene3"; transcript_id "transcript4";'
        #print(result)
        self.assertEqual(str(result).strip(), expected_result)

        bed  = "./tests/test_data/test.bed"
        result = preprocess_file(bed)
        expected_result =   'test_chr	50	100\n'\
                            'test_chr	150	200\n'\
                            'test_chr	200	245'
        self.assertEqual(str(result).strip(), expected_result)


    def test_merge_targets_merge(self):
        '''known issue with bedtools cat https://github.com/daler/pybedtools/issues/282'''

        input = ["./tests/test_data/test.gtf", "./tests/test_data/test.bed"]
        result = merge_targets(files=input, gtf_feature="transcript", operation="merge", window = [0,0])
        expected_result =   'test_chr	19	100\n'\
                            'test_chr	119	280'
        self.assertEqual(str(result).strip(), expected_result)

        result = merge_targets(files=input, gtf_feature="CDS", operation="intersect", window = [0,0])
        expected_result =   'test_chr	150	160\n'\
                            'test_chr	210	230\n'\
                            'test_chr	239	245'
        self.assertEqual(str(result).strip(), expected_result)

        result = merge_targets(files=["./tests/test_data/test.bed"], window = [5,7])
        expected_result =   'test_chr	45	107\n'\
                            'test_chr	145	252'
        
        self.assertEqual(str(result).strip(), expected_result)

        print("\n\n\ttest_merge_targets_merge triggers a known issue in pybedtools, https://github.com/daler/pybedtools/issues/282 \
            \n\toutput is correct and all tests should still pass.\n")


    def test_gtf_to_tss_bed(self):
        input_gtf = "./tests/test_data/test.gtf"  # Replace with the path to your test GTF file
        tss_upstream = 50
        tss_downstream = 25
        tes_upstream = 5
        tes_downstream = 10

        expected_TSS = [
            'test_chr\t0\t45\tTSS_transcript1\t0\t+',
            'test_chr\t69\t145\tTSS_transcript2\t0\t+',
            'test_chr\t234\t310\tTSS_transcript3\t0\t-',
            'test_chr\t254\t330\tTSS_transcript4\t0\t-',
        ]
        expected_TES = [
            'test_chr\t45\t60\tTES_transcript1\t0\t+',
            'test_chr\t165\t180\tTES_transcript2\t0\t+',
            'test_chr\t190\t204\tTES_transcript3\t0\t-',
            'test_chr\t190\t204\tTES_transcript4\t0\t-'
        ]


        result_TSS, result_TES = gtf_to_tss_tes_bed(input_gtf, 
                                                    tss_upstream=tss_upstream, tss_downstream=tss_downstream,
                                                    tes_upstream=tes_upstream, tes_downstream=tes_downstream)

        self.assertEqual(result_TSS, expected_TSS)
        self.assertEqual(result_TES, expected_TES)

class TestParseLineFunction(unittest.TestCase):

    def test_parse_gtf_line(self):
        line = 'chr1\tENSEMBL\tgene\t11869\t14412\t.\t+\t.\tgene_id "ENSG00000223972"; gene_type "transcribed_unprocessed_pseudogene";'
        fields, attributes = parse_line(line)

        expected_fields = ['chr1', 'ENSEMBL', 'gene', '11869', '14412', '.', '+', '.', 'gene_id "ENSG00000223972"; gene_type "transcribed_unprocessed_pseudogene";']
        expected_attributes = {'gene_id': 'ENSG00000223972', 'gene_type': 'transcribed_unprocessed_pseudogene'}

        self.assertEqual(fields, expected_fields)
        self.assertEqual(attributes, expected_attributes)

    def test_parse_gff_line(self):
        line = 'chr1\tENSEMBL\tgene\t11869\t14412\t.\t+\t.\tID=gene:ENSG00000223972; Name=DDX11L1'
        fields, attributes = parse_line(line)

        expected_fields = ['chr1', 'ENSEMBL', 'gene', '11869', '14412', '.', '+', '.', 'ID=gene:ENSG00000223972; Name=DDX11L1']
        expected_attributes = {'ID': 'gene:ENSG00000223972', 'Name': 'DDX11L1'}

        self.assertEqual(fields, expected_fields)
        self.assertEqual(attributes, expected_attributes)

class TestExtractTranscriptGeneRelationship(unittest.TestCase):

    def test_extract_relationship(self):
        input_file = './tests/test_data/test.gtf'
        expected_relationship = {
            'transcript1': 'gene1',
            'transcript2': 'gene2',
            'transcript3': 'gene3',
            'transcript4': 'gene3'
        }

        result = extract_transcript_gene_relationship(input_file)
        self.assertEqual(result, expected_relationship)

class TestParseInput(unittest.TestCase):

    def test_parse_input(self):
        # Mock input GTF data as a string to simulate file reading
        input_file = './tests/test_data/test.gtf'

        # Expected result format based on the provided GTF data
        expected_genes = {
            'gene1': {
                'exons': {
                    'transcript1': [(20, 35), (40, 50)]
                },
                'CDS_coords': {
                    'transcript1': [(25, 30)]
                },
                'CDS': set(['transcript1']),
                'strand': '+',
                'chromosome': 'test_chr'
            },
            'gene2': {
                'exons': {
                    'transcript2': [(120, 140), (150, 170)]
                },
                'CDS_coords': {
                    'transcript2': [(130, 140), (150, 160)]
                },
                'CDS': set(['transcript2']),
                'strand': '+',
                'chromosome': 'test_chr'
            },
            'gene3': {
                'exons': {
                    'transcript3': [(200, 230), (240, 260)],
                    'transcript4': [(200, 230), (240, 260)]
                },
                'CDS_coords': {
                    'transcript3': [(211, 230), (240, 251)],
                    'transcript4': [(211, 230), (240, 251)]
                },
                'CDS': set(['transcript3', 'transcript4']),
                'strand': '-',
                'chromosome': 'test_chr'
            }
        }

        result = parse_input(input_file)

        self.assertEqual(result, expected_genes)
        

class TestCreateConstitutiveModel(unittest.TestCase):

    def test_create_constitutive_model(self):
        # Path to your test GTF file
        input_file = "./tests/test_data/test.gtf"  
        output_file = "./tests/test_output/test.consensus.gtf"
        output_str, genes_without_consensus = create_constitutive_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        # Assertions to verify the correctness of the output
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
        self.assertEqual(genes_without_consensus, set())
    
    def test_create_coding_a2b1_model(self):
        # a2b1 coding transcripts
        input_file = "./tests/test_data/Hnrnpa2b1_coding.gtf"  
        output_file = "./tests/test_output/Hnrnpa2b1_coding.consensus.gtf"
        output_str, genes_without_consensus = create_constitutive_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        # Assertions to verify the correctness of the output
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
        self.assertEqual(genes_without_consensus, set())

    def test_create_noncoding_a2b1_model(self):
        # a2b1 non_coding transcripts
        input_file = "./tests/test_data/Hnrnpa2b1_non_coding.gtf"  
        output_file = "./tests/test_output/Hnrnpa2b1_non_coding.consensus.gtf"
        output_str, genes_without_consensus = create_constitutive_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        # Assertions to verify the correctness of the output
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
        self.assertEqual(genes_without_consensus, set())
    
    def test_create_a2b1_model(self):
        # a2b1 all transcripts
        input_file = "./tests/test_data/Hnrnpa2b1.gtf"  
        output_file = "./tests/test_output/Hnrnpa2b1.consensus.gtf"
        output_str, genes_without_consensus = create_constitutive_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        # Assertions to verify the correctness of the output
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
        self.assertEqual(genes_without_consensus, set())
    
    def test_create_a2b1_against_coding_model(self):
        # a2b1 all transcripts should match coding gtf output
        input_file = "./tests/test_data/Hnrnpa2b1.gtf"  
        output_file = "./tests/test_output/Hnrnpa2b1_coding.consensus.gtf"
        output_str, genes_without_consensus = create_constitutive_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        # Assertions to verify the correctness of the output
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
        self.assertEqual(genes_without_consensus, set())

    def test_create_spi1_model(self):
        # a2b1 all transcripts should match coding gtf output
        input_file = "./tests/test_data/Spi1.gtf"  
        output_file = "./tests/test_output/Spi1.consensus.gtf"
        output_str, genes_without_consensus = create_constitutive_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        # Assertions to verify the correctness of the output
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
        self.assertEqual(genes_without_consensus, set())

class TestCreateMetageneModel(unittest.TestCase):

    def test_create_constitutive_model(self):
        # Path to your test GTF file
        input_file = "./tests/test_data/test.gtf"  
        output_file = "./tests/test_output/test.metagene.gtf"
        output_str = create_metagene_model(input_file)        
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
    
    def test_create_coding_a2b1_model(self):
        # a2b1 coding transcripts
        input_file = "./tests/test_data/Hnrnpa2b1_coding.gtf"  
        output_file = "./tests/test_output/Hnrnpa2b1_coding.metagene.gtf"
        output_str = create_metagene_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
 
    def test_create_coding_a2b1_model(self):
        # a2b1 noncoding transcripts
        input_file = "./tests/test_data/Hnrnpa2b1_non_coding.gtf"  
        output_file = "./tests/test_output/Hnrnpa2b1_non_coding.metagene.gtf"
        output_str = create_metagene_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
 
    def test_create_coding_a2b1_model(self):
        # a2b1 coding transcripts
        input_file = "./tests/test_data/Hnrnpa2b1.gtf"  
        output_file = "./tests/test_output/Hnrnpa2b1.metagene.gtf"
        output_str = create_metagene_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))
 

    def test_create_spi1_model(self):
        input_file = "./tests/test_data/Spi1.gtf"  
        output_file = "./tests/test_output/Spi1.metagene.gtf"
        output_str = create_metagene_model(input_file)
        with open(output_file, 'r') as f: expected_output_str = f.read().strip()
        self.assertEqual(sorted(output_str.split("\n")), sorted(expected_output_str.split("\n")))

def create_expected_exons_dataframe(expected_exons, feature='exon'):
    """
    Convert a list of exon dictionaries to a dataframe with the correct headers.
    
    Args:
        expected_exons (list of dict): List of exons with keys "chrom", "start", "end", "strand", and "transcript_id".
        
    Returns:
        pd.DataFrame: DataFrame with the correct headers.
    """
    # Convert to dataframe
    df = pd.DataFrame(expected_exons)

    # Add missing columns with placeholder values
    df['source'] = 'Test'  # Placeholder value
    df['feature'] = feature  # Placeholder value
    df['score'] = '.'  # Placeholder value
    df['frame'] = '.'  # Placeholder value
    df['attributes'] = df.apply(lambda row: f'gene_id "{row.gene_id}"; transcript_id "{row.transcript_id}";', axis=1)

    # Rearrange columns to match the desired header
    df = df[['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']]
    
    return df

class TestTruncateGTF(unittest.TestCase):
    def setUp(self):
        # Create a sample GTF file content
        
        self.gtf_content_1 = """\
test_chr	Test	gene	20	50	.	+	.	gene_id "gene1";
test_chr	Test	transcript	20	50	.	+	.	gene_id "gene1"; transcript_id "transcript1";
test_chr	Test	exon	20	35	.	+	.	gene_id "gene1"; transcript_id "transcript1";
test_chr	Test	exon	40	50	.	+	.	gene_id "gene1"; transcript_id "transcript1";
test_chr	Test	CDS	25	30	.	+	0	gene_id "gene1"; transcript_id "transcript1";
test_chr	Test	5UTR	20	24	.	+	.	gene_id "gene1"; transcript_id "transcript1";
test_chr	Test	3UTR	31	35	.	+	.	gene_id "gene1"; transcript_id "transcript1";
test_chr	Test	3UTR	40	50	.	+	.	gene_id "gene1"; transcript_id "transcript1";
"""
        self.gtf_content_2 = """\
test_chr	Test	gene	120	170	.	+	.	gene_id "gene2";
test_chr	Test	transcript	120	170	.	+	.	gene_id "gene2"; transcript_id "transcript2";
test_chr	Test	exon	120	140	.	+	.	gene_id "gene2"; transcript_id "transcript2";
test_chr	Test	exon	150	170	.	+	.	gene_id "gene2"; transcript_id "transcript2";
test_chr	Test	CDS	130	140	.	+	.	gene_id "gene2"; transcript_id "transcript2";
test_chr	Test	CDS	150	160	.	+	.	gene_id "gene2"; transcript_id "transcript2";
test_chr	Test	5UTR	120	129	.	+	.	gene_id "gene2"; transcript_id "transcript2";
test_chr	Test	3UTR	161	170	.	+	.	gene_id "gene2"; transcript_id "transcript2";
"""
        self.gtf_content_3 = """\
test_chr	Test	gene	200	280	.	-	.	gene_id "gene3";
test_chr	Test	transcript	200	260	.	-	.	gene_id "gene3"; transcript_id "transcript3";
test_chr	Test	exon	200	230	.	-	.	gene_id "gene3"; transcript_id "transcript3";
test_chr	Test	exon	240	260	.	-	.	gene_id "gene3"; transcript_id "transcript3";
test_chr	Test	CDS	211	230	.	-	0	gene_id "gene3"; transcript_id "transcript3";
test_chr	Test	CDS	240	251	.	-	0	gene_id "gene3"; transcript_id "transcript3";
test_chr	Test	5UTR	252	260	.	-	.	gene_id "gene3"; transcript_id "transcript3";
test_chr	Test	3UTR	200	210	.	-	.	gene_id "gene3"; transcript_id "transcript3";
test_chr	Test	transcript	200	280	.	-	.	gene_id "gene3"; transcript_id "transcript4";
test_chr	Test	exon	200	230	.	-	.	gene_id "gene3"; transcript_id "transcript4";
test_chr	Test	exon	240	260	.	-	.	gene_id "gene3"; transcript_id "transcript4";
test_chr	Test	CDS	211	230	.	-	0	gene_id "gene3"; transcript_id "transcript4";
test_chr	Test	CDS	240	251	.	-	0	gene_id "gene3"; transcript_id "transcript4";
test_chr	Test	5UTR	252	280	.	-	.	gene_id "gene3"; transcript_id "transcript4";
test_chr	Test	3UTR	200	210	.	-	.	gene_id "gene3"; transcript_id "transcript4";        
"""

        self.gtf_content_4 = """\
test_chr	Test	gene	200	330	.	-	.	gene_id "gene3";
test_chr	Test	transcript	200	330	.	-	.	gene_id "gene3"; transcript_id "transcript3";
test_chr	Test	exon	200	230	.	-	.	gene_id "gene3"; transcript_id "transcript3";
test_chr	Test	exon	240	260	.	-	.	gene_id "gene3"; transcript_id "transcript3"; 
test_chr	Test	exon	280	330	.	-	.	gene_id "gene3"; transcript_id "transcript3";    
"""
        self.gtf_df = pd.read_csv(StringIO(self.gtf_content_1), sep='\t', header=None, comment='#')
        self.gtf_df.columns = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]

   
    def test_truncate_gtf_no_truncation(self):
        df, unique_transcript_ids, unique_gene_ids = truncate_gtf(StringIO(self.gtf_content_1), feature="exon", percentiles=[0,100])
        self.assertEqual(len(df), len(self.gtf_df[self.gtf_df["feature"] == "exon"]))

        expected_exons = [
            {"chrom": "test_chr", "start": 20, "end": 35, "strand": "+", "gene_id": "gene1", "transcript_id": "transcript1"},
            {"chrom": "test_chr", "start": 40, "end": 50, "strand": "+", "gene_id": "gene1", "transcript_id": "transcript1"},
        ]

        expected_df = create_expected_exons_dataframe(expected_exons)

        expected_dict = expected_df.to_dict(orient='list')
        test_dict = df.to_dict(orient='list')

        self.assertEqual(expected_dict, test_dict)

        df, unique_transcript_ids, unique_gene_ids = truncate_gtf(StringIO(self.gtf_content_1), feature="exon", percentiles=[0,50])

        expected_exons = [
            {"chrom": "test_chr", "start": 20, "end": 34, "strand": "+", "gene_id": "gene1", "transcript_id": "transcript1"},
        ]

        expected_df = create_expected_exons_dataframe(expected_exons)
        expected_dict = expected_df.to_dict(orient='list')
        test_dict = df.to_dict(orient='list')

        self.assertEqual(expected_dict, test_dict)

        df, unique_transcript_ids, unique_gene_ids = truncate_gtf(StringIO(self.gtf_content_1), feature="exon", percentiles=[10,20])

        expected_exons = [
            {"chrom": "test_chr", "start": 23, "end": 25, "strand": "+", "gene_id": "gene1", "transcript_id": "transcript1"},
        ]

        expected_df = create_expected_exons_dataframe(expected_exons)
        expected_dict = expected_df.to_dict(orient='list')
        test_dict = df.to_dict(orient='list')

        self.assertEqual(expected_dict, test_dict)

        ##############################
        ############################## gtf_2
        ##############################

        df, unique_transcript_ids, unique_gene_ids = truncate_gtf(StringIO(self.gtf_content_2), feature="CDS", percentiles=[25,75])

        expected_exons = [
            {"chrom": "test_chr", "start": 136, "end": 140, "strand": "+", "gene_id": "gene2", "transcript_id": "transcript2"},
            {"chrom": "test_chr", "start": 150, "end": 155, "strand": "+", "gene_id": "gene2", "transcript_id": "transcript2"},
        ]

        expected_df = create_expected_exons_dataframe(expected_exons, feature="CDS")
        expected_dict = expected_df.to_dict(orient='list')
        test_dict = df.to_dict(orient='list')

        self.assertEqual(expected_dict, test_dict)

        ##############################
        ############################## gtf_3
        ##############################


        df, unique_transcript_ids, unique_gene_ids = truncate_gtf(StringIO(self.gtf_content_3), feature="exon", percentiles=[0,100])

        expected_exons = [
            {"chrom": "test_chr", "start": 200, "end": 230, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript3"},
            {"chrom": "test_chr", "start": 200, "end": 230, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript4"},
            {"chrom": "test_chr", "start": 240, "end": 260, "strand": "-", "gene_id":  "gene3", "transcript_id": "transcript3"}, 
            {"chrom": "test_chr", "start": 240, "end": 260, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript4"}
        ]

        expected_df = create_expected_exons_dataframe(expected_exons)
        expected_df = expected_df.sort_values(["chrom", "start"])
        expected_dict = expected_df.to_dict(orient='list')
        test_dict = df.to_dict(orient='list')

        self.assertEqual(expected_dict, test_dict)


        df, unique_transcript_ids, unique_gene_ids = truncate_gtf(StringIO(self.gtf_content_3), feature="exon", percentiles=[0,50])


        expected_exons = [
            {"chrom": "test_chr", "start": 226, "end": 230, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript3"},
            {"chrom": "test_chr", "start": 226, "end": 230, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript4"},
            {"chrom": "test_chr", "start": 240, "end": 260, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript3"}, 
            {"chrom": "test_chr", "start": 240, "end": 260, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript4"}
        ]
        expected_df = create_expected_exons_dataframe(expected_exons)
        expected_df = expected_df.sort_values(["chrom", "start"])
        expected_dict = expected_df.to_dict(orient='list')
        test_dict = df.to_dict(orient='list')

        self.assertEqual(expected_dict, test_dict)

        df, unique_transcript_ids, unique_gene_ids = truncate_gtf(StringIO(self.gtf_content_4), feature="exon", percentiles=[25,75])


        expected_exons = [
            {"chrom": "test_chr", "start": 226, "end": 230, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript3"},
            {"chrom": "test_chr", "start": 240, "end": 260, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript3"},
            {"chrom": "test_chr", "start": 280, "end": 305, "strand": "-", "gene_id": "gene3", "transcript_id": "transcript3"}

        ]

        expected_df = create_expected_exons_dataframe(expected_exons)
        expected_df = expected_df.sort_values(["chrom", "start"])
        expected_dict = expected_df.to_dict(orient='list')
        test_dict = df.to_dict(orient='list')

        self.assertEqual(expected_dict, test_dict)

if __name__ == '__main__':
    unittest.main()