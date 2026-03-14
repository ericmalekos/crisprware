import unittest
import subprocess
import sys


class TestCLI(unittest.TestCase):
    def test_version_flag(self):
        result = subprocess.run(
            [sys.executable, "-m", "crisprware.cli", "-v"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("crisprware", result.stdout)

    def test_no_args_exits_nonzero(self):
        result = subprocess.run(
            [sys.executable, "-m", "crisprware.cli"],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_help_lists_subcommands(self):
        result = subprocess.run(
            [sys.executable, "-m", "crisprware.cli", "--help"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("generate_guides", result.stdout)
        self.assertIn("score_guides", result.stdout)
        self.assertIn("rank_guides", result.stdout)


if __name__ == "__main__":
    unittest.main()
