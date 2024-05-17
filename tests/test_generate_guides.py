import unittest
from argparse import Namespace
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from unittest.mock import patch
from src.generate_guides import parse_arguments, find_sgRNA, output_bed_line, process_pam, reverse_cut_site_offset
from utils.dna_sequence_functions import revcom

class TestGenerateGuides(unittest.TestCase):

    # @patch('sys.argv', ['generate_guides.py', '-f', 'test.fa', '-p', 'NGG', '-l', '20'])
    # def test_parse_arguments(self):
    #     args = parse_arguments()
    #     self.assertEqual(args.fasta, "test.fa")
    #     self.assertEqual(args.pam, "NGG")
    #     self.assertEqual(args.sgRNA_length, 20)

    def test_find_sgRNA(self):
        args = Namespace(context_window=(4, 6), sgRNA_length=20, pam_5_prime=False)
        pam = "CGG"
        chrm = """CCGTACGATCCGACGTTTGACCTGATCCGTACGGATCCGGGATGCC"""
        start = 0
        end = len(chrm)
        
        sgRNAs = list(find_sgRNA(args, pam, chrm, start, end))
                
        expected_sgRNAs = [
            ('GACGTTTGACCTGATCCGTA', 12, 'ATCCGACGTTTGACCTGATCCGTACGGATC'),
            ('TGACCTGATCCGTACGGATC', 18, 'CGTTTGACCTGATCCGTACGGATCCGGGAT')
        ]
        self.assertEqual(sgRNAs, expected_sgRNAs)

        pam = "CCT"
        sgRNAs = list(find_sgRNA(args, pam, chrm, start, end, forward=False))

        expected_sgRNAs = [('GATCCGTACGGATCCGGGAT', 21, 'TGACCTGATCCGTACGGATCCGGGATGCC')]

        self.assertEqual(sgRNAs, expected_sgRNAs)

        chrm = """CCGTACGATCCGACGTTTGACCTGATCCGTACGAAAGATCCGGGATGCC"""
        args = Namespace(context_window=(7, 4), sgRNA_length=23, pam_5_prime=True)
        pam = "TTTG"
        expected_sgRNAs = [('ACCTGATCCGTACGAAAGATCCG', 16, 'ACGTTTGACCTGATCCGTACGAAAGATCCGGGAT')]
        sgRNAs = list(find_sgRNA(args, pam, chrm, start, end, forward=True))
        self.assertEqual(sgRNAs, expected_sgRNAs)

        pam = "GAAA"
        sgRNAs = list(find_sgRNA(args, pam, chrm, start, end, forward=False))
        expected_sgRNAs = [('CCGACGTTTGACCTGATCCGTAC', 33, 'CGATCCGACGTTTGACCTGATCCGTACGAAAGAT')]
        self.assertEqual(sgRNAs, expected_sgRNAs)


    def test_output_bed_line(self):
        args = Namespace(prefix="", active_site_offset_5=-4, active_site_offset_3=-4, pam="NGG", pam_5_prime=False)
        chrm_name = "chr1"
        sgRNA = {
            "sequence": "GACGTTTGACCTGATCCGTA",
            "position": 12,
            "pam": "NGG",
            "sense": "+",
            "length": 20,
            "context": "ATCCGACGTTTGACCTGATCCGTACGGATC"
        }

        bed_line = output_bed_line(args, chrm_name, sgRNA)
        expected_bed_line = "chr1\t27\t27\tchr1:12:+,GACGTTTGACCTGATCCGTA,NGG,chr1,12,+\tATCCGACGTTTGACCTGATCCGTACGGATC\t+"
        self.assertEqual(bed_line, expected_bed_line)


        sequence = revcom("GATCCGTACGGATCCGGGAT")
        context = revcom("TGACCTGATCCGTACGGATCCGGGATGCC")
        sgRNA = {
            "sequence": sequence,
            "position": 23,
            "pam": "NGG",
            "sense": "-",
            "length": 20,
            "context": context
        }

        bed_line = output_bed_line(args, chrm_name, sgRNA)
        expected_bed_line = "chr1\t29\t29\tchr1:23:-,"+sequence+",NGG,chr1,23,-\t"+context+"\t-"
        self.assertEqual(bed_line, expected_bed_line)

        sequence = "ACCTGATCCGTACGAAAGATCCG"
        context = "ACGTTTGACCTGATCCGTACGAAAGATCCGGGAT"
        args = Namespace(prefix="", active_site_offset_5=19, active_site_offset_3=23, pam="TTTV", pam_5_prime=True)
        sgRNA = {
            "sequence": "ACCTGATCCGTACGAAAGATCCG",
            "position": 23,
            "pam": "TTTV",
            "sense": "+",
            "length": 23,
            "context": "ACGTTTGACCTGATCCGTACGAAAGATCCGGGAT"
        }
        bed_line = output_bed_line(args, chrm_name, sgRNA)
        expected_bed_line = "chr1\t45\t49\tchr1:23:+,"+sequence+",TTTV,chr1,23,+\t"+context+"\t+"
        self.assertEqual(bed_line, expected_bed_line)

        sequence = revcom("CCGACGTTTGACCTGATCCGTAC")
        context = revcom("CGATCCGACGTTTGACCTGATCCGTACGAAAGAT")
        args = Namespace(prefix="", active_site_offset_5=19, active_site_offset_3=23, pam="TTTV", pam_5_prime=True)
        sgRNA = {
            "sequence": sequence,
            "position": 33,
            "pam": "TTTV",
            "sense": "-",
            "length": 23,
            "context": context
        }
        bed_line = output_bed_line(args, chrm_name, sgRNA)
        expected_bed_line = "chr1\t9\t13\tchr1:33:-,"+sequence+",TTTV,chr1,33,-\t"+context+"\t-"
        self.assertEqual(bed_line, expected_bed_line)


    def test_reverse_cut_site_offset(self):
        
        args = Namespace(prefix="", active_site_offset_5=-4, active_site_offset_3=-4, sgRNA_length=20, pam_5_prime=False, pam="AGG")
        bed_line = "chr10\t7163868\t7163868\tchr10:7163852:+,TCCATGCACTCCTGCACCGT,NGG,chr10,7163852,+\tATCCGACGTTTGACCTGATCCGTACGGATC\t+"
        bed_line = reverse_cut_site_offset(bed_line, args)
        expected_bed_line="chr10\t7163851\t7163871\tchr10:7163852:+,TCCATGCACTCCTGCACCGT,NGG,chr10,7163852,+\tATCCGACGTTTGACCTGATCCGTACGGATC\t+"
        self.assertEqual(bed_line, expected_bed_line)


        bed_line = "chr10\t7163802\t7163802\tchr10:7163796:-,CCACTTTGACTCTCTCGATC,NGG,chr10,7163796,-\tATCCGACGTTTGACCTGATCCGTACGGATC\t-"
        bed_line = reverse_cut_site_offset(bed_line, args)
        expected_bed_line="chr10\t7163798\t7163818\tchr10:7163796:-,CCACTTTGACTCTCTCGATC,NGG,chr10,7163796,-\tATCCGACGTTTGACCTGATCCGTACGGATC\t-"
        self.assertEqual(bed_line, expected_bed_line)


        args = Namespace(prefix="", active_site_offset_5=19, active_site_offset_3=23, sgRNA_length=23, pam="TTTV", pam_5_prime=True)

        bed_line = "chr10\t7163868\t7163868\tchr10:7163852:+,TCCATGCACTCCTGCACCGT,TTTV,chr10,7163852,+\tATCCGACGTTTGACCTGATCCGTACGGATC\t+"
        bed_line = reverse_cut_site_offset(bed_line, args)
        expected_bed_line="chr10\t7163855\t7163878\tchr10:7163852:+,TCCATGCACTCCTGCACCGT,TTTV,chr10,7163852,+\tATCCGACGTTTGACCTGATCCGTACGGATC\t+"
        self.assertEqual(bed_line, expected_bed_line) 

        bed_line = "chr10\t7163868\t7163868\tchr10:7163852:-,TCCATGCACTCCTGCACCGT,TTTV,chr10,7163852,-\tATCCGACGTTTGACCTGATCCGTACGGATC\t-"
        bed_line = reverse_cut_site_offset(bed_line, args)
        expected_bed_line="chr10\t7163828\t7163851\tchr10:7163852:-,TCCATGCACTCCTGCACCGT,TTTV,chr10,7163852,-\tATCCGACGTTTGACCTGATCCGTACGGATC\t-"
        self.assertEqual(bed_line, expected_bed_line) 

    # def test_process_pam(self):
    #     args = Namespace(context_window=(4, 6), sgRNA_length=20, pam_5_prime=False)
    #     pam = "NGG"
    #     record = SeqRecord(Seq("ACGTACGTNGGACGTACGTNGGACGTACGT"), id="chr1")
    #     start = 0
    #     end = len(record.seq)
    #     pam_set = {"NGG"}
    #     rev_pam_set = {"CCN"}

    #     results = process_pam(args, pam, record, start, end, pam_set, rev_pam_set)
    #     expected_results = [
    #         "chr1\t16\t18\tchr1:6:+,ACGTACGTNGGACGTACGTNG,NGG,chr1,6,+\tACGTACGTNGGACGTACGTNGGACGTACGT\t+\n",
    #         "chr1\t27\t29\tchr1:17:+,ACGTACGTNGGACGTACGTNG,NGG,chr1,17,+\tACGTACGTNGGACGTACGTNGGACGTACGT\t+\n"
    #     ]
    #     self.assertEqual(results, expected_results)

if __name__ == "__main__":
    unittest.main()