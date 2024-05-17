import unittest
from utils.gtf_bed_processing_functions import preprocess_file,merge_targets,\
    gtf_to_tss_tes_bed,parse_line, extract_transcript_gene_relationship, parse_input,\
        create_metagene_model,create_constitutive_model, truncate_gtf
import pandas as pd

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
 

# class TestTruncateGTF(unittest.TestCase):
#     def setUp(self):
#         # This method is called before each test
#         self.input_file = 'tests/test_data/test.gtf'
#         self.feature = 'CDS'
#         self.percentiles = [0, 50]

#     def test_equivalent_output(self):
#         # Process the test GTF with both functions
#         df1 = truncate_gtf(self.input_file, self.feature, self.percentiles)

#         df2 = truncate_gtf_2(self.input_file, self.feature, self.percentiles)

#         print("\t\ttruncate_gtf")
#         print(df1)

#         print("\t\ttruncate_gtf_2")
#         print(df2)

#         # Verify the DataFrames are equivalent
#         pd.testing.assert_frame_equal(df1, df2, check_dtype=True, check_like=True)

if __name__ == '__main__':
    unittest.main()