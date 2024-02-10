import unittest
import pandas as pd
from utils.quantified_rna_functions import calculate_statistics,process_kallisto,\
process_salmon,process_flair,process_mandalorian,process_dataframes,infer_file_type_from_first_line,\
process_files,filter_dataframe,add_gene_ids_and_subset


class TestProcessRNA(unittest.TestCase):

    def test_process_flair_1(self):
        flair_counts = ["./tests/test_data/processed_rna_seq/flair/flair1_counts_matrix.tsv"]
        output = process_flair(flair_counts)
        expected_output = pd.read_csv("./tests/test_output/processed_rna_out/flair_1_output.tsv", delimiter='\t')
        self.assertTrue(output.equals(expected_output), "DataFrames are not equal")


    def test_process_flair_2(self):
        flair_counts = ["./tests/test_data/processed_rna_seq/flair/flair2_counts_matrix.tsv",\
                        "./tests/test_data/processed_rna_seq/flair/flair2_counts_matrix.tsv"]
        output = process_flair(flair_counts)
        expected_output = pd.read_csv("./tests/test_output/processed_rna_out/flair_2_output.tsv", delimiter='\t')
        self.assertTrue(output.equals(expected_output), "DataFrames are not equal")

###############################################################################################################
        
    def test_process_kallisto_1(self):
        kallisto_tpms = ["./tests/test_data/processed_rna_seq/kallisto/kallisto_abundance_1.tsv"]
        output = process_kallisto(kallisto_tpms)
        expected_output = pd.read_csv("./tests/test_output/processed_rna_out/kallisto_1_output.tsv", delimiter='\t')
        self.assertTrue(output.equals(expected_output), "DataFrames are not equal")

    def test_process_kallisto_2(self):
        kallisto_tpms = ["./tests/test_data/processed_rna_seq/kallisto/kallisto_abundance_1.tsv",
                        "./tests/test_data/processed_rna_seq/kallisto/kallisto_abundance_2.tsv"]
        output = process_kallisto(kallisto_tpms)
        expected_output = pd.read_csv("./tests/test_output/processed_rna_out/kallisto_2_output.tsv", delimiter='\t')
        self.assertTrue(output.equals(expected_output), "DataFrames are not equal")

###############################################################################################################

    def test_process_salmon_1(self):
        salmon_tpms = ["./tests/test_data/processed_rna_seq/salmon/salmon_quant_1.sf"]
        output = process_salmon(salmon_tpms)
        expected_output = pd.read_csv("./tests/test_output/processed_rna_out/salmon_1_output.tsv", delimiter='\t')
        self.assertTrue(output.equals(expected_output), "DataFrames are not equal")

    def test_process_salmon_2(self):
        salmon_tpms = ["./tests/test_data/processed_rna_seq/salmon/salmon_quant_1.sf",
                       "./tests/test_data/processed_rna_seq/salmon/salmon_quant_2.sf"]
        output = process_salmon(salmon_tpms)
        #output.to_csv("./tests/test_output/processed_rna_out/salmon_2_output.tsv", sep="\t", index=False)
        expected_output = pd.read_csv("./tests/test_output/processed_rna_out/salmon_2_output.tsv", delimiter='\t')
        self.assertTrue(output.equals(expected_output), "DataFrames are not equal")

###############################################################################################################
    
    def test_process_mandalorian_1(self):
        mandalorian_tpms = ["./tests/test_data/processed_rna_seq/mandalorian/mandalorian_isoforms_1.filtered.clean.tpm"]
        output = process_mandalorian(mandalorian_tpms)
        output.to_csv("./tests/test_output/processed_rna_out/mandalorian_1_output.tsv", sep="\t", index=False)
        expected_output = pd.read_csv("./tests/test_output/processed_rna_out/mandalorian_1_output.tsv", delimiter='\t')
        self.assertTrue(output.equals(expected_output), "DataFrames are not equal")

    def test_process_mandalorian_2(self):
        mandalorian_tpms = ["./tests/test_data/processed_rna_seq/mandalorian/mandalorian_isoforms_1.filtered.clean.tpm",
                       "./tests/test_data/processed_rna_seq/mandalorian/mandalorian_isoforms_2.filtered.clean.tpm"]
        output = process_mandalorian(mandalorian_tpms)
        expected_output = pd.read_csv("./tests/test_output/processed_rna_out/mandalorian_2_output.tsv", delimiter='\t')
        self.assertTrue(output.equals(expected_output), "DataFrames are not equal")

if __name__ == '__main__':
    unittest.main()