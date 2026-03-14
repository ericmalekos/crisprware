import unittest
import tempfile
import os
import shutil

import pandas as pd

from crisprware.utils.rank_guides_functions import (
    create_combined_weighted_column,
    filter_df_by_column,
    group_and_minimize,
    select_guides,
    validate_and_modify_bed,
    analyze_target_ids,
)


class TestCreateCombinedWeightedColumn(unittest.TestCase):
    def test_equal_weights(self):
        df = pd.DataFrame({"score_a": [0.0, 0.5, 1.0], "score_b": [1.0, 0.5, 0.0]})
        result = create_combined_weighted_column(df, ["score_a", "score_b"])
        self.assertIn("combined_weighted", result.columns)
        self.assertIn("score_a_normalized", result.columns)
        self.assertIn("score_b_normalized", result.columns)

    def test_custom_weights(self):
        df = pd.DataFrame({"a": [0.0, 10.0], "b": [0.0, 10.0]})
        result = create_combined_weighted_column(df, ["a", "b"], weights=[2.0, 0.0])
        # Weight of 0 on b means only a matters
        self.assertEqual(result["combined_weighted"].iloc[0], 0.0)

    def test_mismatched_lengths(self):
        df = pd.DataFrame({"a": [1], "b": [2]})
        with self.assertRaises(ValueError):
            create_combined_weighted_column(df, ["a", "b"], weights=[1.0])


class TestFilterDfByColumn(unittest.TestCase):
    def test_filters_below_cutoff(self):
        df = pd.DataFrame({"val": [1, 5, 10, 15]})
        result = filter_df_by_column(df, "val", 6)
        self.assertEqual(len(result), 2)
        self.assertTrue((result["val"] >= 6).all())


class TestGroupAndMinimize(unittest.TestCase):
    def test_keep_top_n(self):
        df = pd.DataFrame(
            {
                "target_id": ["A", "A", "A", "B", "B"],
                "rank": [3, 1, 2, 5, 4],
            }
        )
        result = group_and_minimize(df, "rank", 2)
        # Should keep 2 per target
        a_count = (result["target_id"] == "A").sum()
        b_count = (result["target_id"] == "B").sum()
        self.assertEqual(a_count, 2)
        self.assertEqual(b_count, 2)

    def test_keep_all(self):
        df = pd.DataFrame({"target_id": ["A", "A", "B"], "rank": [1, 2, 3]})
        result = group_and_minimize(df, "rank", -1)
        self.assertEqual(len(result), 3)


class TestSelectGuides(unittest.TestCase):
    def test_no_spacing(self):
        df = pd.DataFrame(
            {
                "target_id": ["A", "A"],
                "start": [100, 105],
                "stop": [120, 125],
                "rank": [2, 1],
            }
        )
        result = select_guides(df, "rank", 0)
        self.assertEqual(len(result), 2)

    def test_with_spacing(self):
        df = pd.DataFrame(
            {
                "target_id": ["A", "A", "A"],
                "start": [100, 105, 200],
                "stop": [120, 125, 220],
                "rank": [3, 2, 1],
            }
        )
        result = select_guides(df, "rank", 50)
        # First (rank=3, 100-120) kept, second (105-125) too close, third (200-220) kept
        self.assertEqual(len(result), 2)


class TestValidateAndModifyBed(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_valid_bed(self):
        path = os.path.join(self.tmpdir, "test.bed")
        with open(path, "w") as f:
            f.write("chr1\t100\t200\tgene1\n")
            f.write("chr1\t300\t400\tgene2\n")
        df, target_ids = validate_and_modify_bed(path)
        self.assertEqual(len(df), 2)
        self.assertEqual(target_ids, {"gene1", "gene2"})

    def test_too_few_columns(self):
        path = os.path.join(self.tmpdir, "bad.bed")
        with open(path, "w") as f:
            f.write("chr1\t100\n")
        with self.assertRaises(ValueError):
            validate_and_modify_bed(path)

    def test_start_greater_than_end(self):
        path = os.path.join(self.tmpdir, "bad.bed")
        with open(path, "w") as f:
            f.write("chr1\t500\t100\tgene1\n")
        with self.assertRaises(ValueError):
            validate_and_modify_bed(path)


class TestAnalyzeTargetIds(unittest.TestCase):
    def test_returns_counts(self):
        df = pd.DataFrame({"target_id": ["A", "A", "A", "B", "B"]})
        counts = analyze_target_ids(df, {"C"})
        self.assertEqual(counts["A"], 3)
        self.assertEqual(counts["B"], 2)


if __name__ == "__main__":
    unittest.main()
