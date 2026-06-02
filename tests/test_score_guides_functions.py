import unittest
import argparse
import tempfile
import os
import shutil

from crisprware.score_guides import restricted_float, check_files_exist, get_alt_pams


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


class TestCheckFilesExist(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_missing_files(self):
        with self.assertRaises(FileNotFoundError):
            check_files_exist(os.path.join(self.tmpdir, "nonexistent"))

    def test_existing_files(self):
        index = os.path.join(self.tmpdir, "genome.fa.index")
        for ext in [".reverse", ".forward", ".gs"]:
            with open(index + ext, "w") as f:
                f.write("")
        # Should not raise
        check_files_exist(index)


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
