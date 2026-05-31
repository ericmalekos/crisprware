import unittest
import argparse
import tempfile
import os
import shutil

from crisprware.score_guides import restricted_float, detect_indexer, gscan_scoring, get_alt_pams


class TestRestrictedFloat(unittest.TestCase):
    def test_valid_values(self):
        self.assertEqual(restricted_float(0.0), 0.0)
        self.assertEqual(restricted_float(0.5), 0.5)
        self.assertEqual(restricted_float(1.0), 1.0)
        self.assertEqual(restricted_float("0.75"), 0.75)

    def test_out_of_range_high(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            restricted_float(1.5)

    def test_out_of_range_low(self):
        with self.assertRaises(argparse.ArgumentTypeError):
            restricted_float(-0.1)


class TestDetectIndexer(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_missing_files(self):
        with self.assertRaises(FileNotFoundError):
            detect_indexer(os.path.join(self.tmpdir, "nonexistent"))

    def test_detects_guidescan2(self):
        index = os.path.join(self.tmpdir, "genome_guidescan2")
        for ext in [".reverse", ".forward", ".gs"]:
            with open(index + ext, "w") as f:
                f.write("")
        self.assertEqual(detect_indexer(index), "guidescan2")

    def test_detects_crispr_ots(self):
        # crispr-ots writes a real .crot alongside zero-byte guidescan-compatible
        # stubs, so the .crot file must take precedence over the stub trio.
        index = os.path.join(self.tmpdir, "genome_crisprots")
        for ext in [".crot", ".reverse", ".forward", ".gs"]:
            with open(index + ext, "w") as f:
                f.write("")
        self.assertEqual(detect_indexer(index), "crispr-ots")


class TestGscanScoringBulgeGuard(unittest.TestCase):
    def test_crispr_ots_rejects_rna_bulges(self):
        # The guard must fire before any subprocess is launched.
        with self.assertRaises(ValueError):
            gscan_scoring(guideCSV="in.csv", output="out.csv", guideIndex="idx", indexer="crispr-ots", rna_bulges=1)

    def test_crispr_ots_rejects_dna_bulges(self):
        with self.assertRaises(ValueError):
            gscan_scoring(guideCSV="in.csv", output="out.csv", guideIndex="idx", indexer="crispr-ots", dna_bulges=2)


class TestGetAltPams(unittest.TestCase):
    def test_single_unambiguous(self):
        result = get_alt_pams(["AGG"])
        self.assertEqual(result, "AGG")

    def test_ambiguous_n(self):
        result = get_alt_pams(["NAG"])
        pams = set(result.split())
        self.assertEqual(pams, {"AAG", "TAG", "CAG", "GAG"})

    def test_multiple_pams(self):
        result = get_alt_pams(["AGG", "NAG"])
        pams = set(result.split())
        self.assertIn("AGG", pams)
        self.assertIn("AAG", pams)
        self.assertIn("TAG", pams)


if __name__ == "__main__":
    unittest.main()
