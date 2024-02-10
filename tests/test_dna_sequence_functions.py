import unittest
from itertools import product
import os
from utils.dna_sequence_functions import map_ambiguous_sequence,subset_fasta_with_bed


class TestIUPAC(unittest.TestCase):
    
    def test_map_ambiguous_sequence(self):
        # Define a test case
        test_sequence = "RYSW"
        # Expected result
        expected_sequences = [
            ''.join(p) for p in product(
                "AG",  # R: AG
                "CT",  # Y: CT
                "GC",  # S: GC
                "AT"   # W: AT
            )
        ]
        
        result = map_ambiguous_sequence(test_sequence)
        self.assertEqual(sorted(result), sorted(expected_sequences))

class TestSubsetFastaWithBed(unittest.TestCase):
    
    def setUp(self):
        # Create sample FASTA and BED files for testing
        self.test_fasta_path = "test.fasta"
        self.test_bed_path = "test.bed"
        self.output_fasta_path = "output.fasta"

        # Create a more varied sample FASTA file content
        with open(self.test_fasta_path, 'w') as file:
            file.write(">chr1\n" + "ATGC" * 25 + "\n")  # 100 bases (25 repeats of ATGC)
            file.write(">chr2\n" + "CGTA" * 25 + "\n")  # Another sequence with different pattern

        # Create sample BED file content
        with open(self.test_bed_path, 'w') as file:
            file.write("chr1\t10\t20\n")  # Valid interval within chr1
            file.write("chr1\t90\t110\n")  # Interval extends beyond sequence length of chr1
            file.write("chr3\t10\t20\n")  # Chromosome not in FASTA

    def test_subset_fasta_with_valid_bed(self):
        # Test subsetting with valid BED intervals
        subset_fasta_with_bed(self.test_fasta_path, self.test_bed_path, self.output_fasta_path)
        with open(self.output_fasta_path) as file:
            content = file.read()
            self.assertIn(">chr1_10:20", content)
            self.assertIn("GCATGCATGC", content)  # Expected sequence from chr1

    def test_bed_interval_extends_beyond_fasta(self):
        # Test handling of BED interval extending beyond FASTA sequence
        subset_fasta_with_bed(self.test_fasta_path, self.test_bed_path, self.output_fasta_path)
        with open(self.output_fasta_path) as file:
            content = file.read()
            self.assertIn(">chr1_90:100", content)
            self.assertIn("GCATGCATGC", content)  # Expected sequence from the end of chr1

    def test_bed_chromosome_not_in_fasta(self):
        # Test handling of BED chromosome not in FASTA
        subset_fasta_with_bed(self.test_fasta_path, self.test_bed_path, self.output_fasta_path)
        with open(self.output_fasta_path) as file:
            content = file.read()
            self.assertNotIn("chr3", content)

    def tearDown(self):
        # Clean up test files
        os.remove(self.test_fasta_path)
        os.remove(self.test_bed_path)
        os.remove(self.output_fasta_path)

if __name__ == '__main__':
    unittest.main()