import unittest
import os
import gzip
import tempfile
import shutil

from crisprware.utils.utility_functions import create_output, decompress_gzip_if_needed, remove_file


class TestCreateOutput(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_with_extension(self):
        full_path, tmp_dir = create_output("/some/path/sample.bed", outdir=self.tmpdir, extension="scored")
        self.assertTrue(full_path.endswith("sample_scored"))
        self.assertIn("sample_scored", os.path.dirname(full_path))
        self.assertTrue(os.path.isdir(os.path.dirname(full_path)))
        self.assertEqual(tmp_dir, "")

    def test_without_extension(self):
        full_path, tmp_dir = create_output("/some/path/sample.bed", outdir=self.tmpdir)
        self.assertTrue(full_path.endswith("sample"))
        self.assertEqual(tmp_dir, "")

    def test_with_stripped(self):
        full_path, _ = create_output("/some/path/sample_gRNA.bed", outdir=self.tmpdir, stripped="_gRNA")
        self.assertIn("sample", full_path)
        self.assertNotIn("_gRNA", full_path)

    def test_with_tmp(self):
        full_path, tmp_dir = create_output("/some/path/sample.bed", outdir=self.tmpdir, extension="out", tmp=True)
        self.assertNotEqual(tmp_dir, "")
        self.assertTrue(os.path.isdir(tmp_dir))
        self.assertTrue(tmp_dir.endswith("tmp/"))


class TestDecompressGzip(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_gz_file(self):
        content = b"hello world"
        gz_path = os.path.join(self.tmpdir, "test.txt.gz")
        with gzip.open(gz_path, "wb") as f:
            f.write(content)

        result_path, was_gzipped = decompress_gzip_if_needed(gz_path)
        self.assertTrue(was_gzipped)
        self.assertEqual(result_path, gz_path[:-3])
        with open(result_path, "rb") as f:
            self.assertEqual(f.read(), content)

    def test_plain_file(self):
        plain_path = os.path.join(self.tmpdir, "test.txt")
        with open(plain_path, "w") as f:
            f.write("hello")

        result_path, was_gzipped = decompress_gzip_if_needed(plain_path)
        self.assertFalse(was_gzipped)
        self.assertEqual(result_path, plain_path)


class TestRemoveFile(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_existing_file(self):
        path = os.path.join(self.tmpdir, "deleteme.txt")
        with open(path, "w") as f:
            f.write("tmp")
        remove_file(path)
        self.assertFalse(os.path.exists(path))

    def test_missing_file(self):
        # Should not raise
        remove_file(os.path.join(self.tmpdir, "nonexistent.txt"))


if __name__ == "__main__":
    unittest.main()
