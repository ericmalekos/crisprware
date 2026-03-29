from __future__ import annotations
import unittest
import os
import tempfile
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from parasol_scripts.score_flashfry_cfd import (
    _split_zero_mismatch_cell,
    _has_ttTV,
    _has_tttt,
    _expected_canonical_coord_from_contig,
    _remove_canonical_zero_mismatch,
    count_mismatches,
    bucket_offtargets_by_mismatch,
    extract_loci_from_cell,
    parse_offtargets,
    build_score_lookup,
    dna_to_rna_base,
    cfd_score_pair,
    format_scores,
    parse_score_list,
    _round_float_str,
    _round_csv_list,
    load_genome_biopython,
    fetch_27mer,
    expand_target_from_contig,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

CHRIII_FASTA = "./tests/test_data/ce11/chrIII_sequence.fasta"
MATRIX_CSV = "./parasol_scripts/off_targ_enCas12a.csv"
FLASHFRY_INPUT = "./tests/test_data/flashfry/ce11_full.part_001.fa.output.tsv"
FLASHFRY_REDUCED = "./tests/test_data/flashfry/ce11_full.part_001.fa.output.reduced.tsv"
UNIQUENESS_INPUT = "./tests/test_data/flashfry/chrIII_uniqueness_test.fa.output.tsv"


def _write_tiny_fasta(path, records: dict[str, str]):
    """Write a dict of {chrom: sequence} to a FASTA file."""
    with open(path, "w") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


def _make_genome_dict(records: dict[str, str]):
    """Return a Biopython SeqIO dict-like from {chrom: seq}."""
    return {name: SeqRecord(Seq(seq), id=name, name=name, description="") for name, seq in records.items()}


# ---------------------------------------------------------------------------
# Pure function tests
# ---------------------------------------------------------------------------


class TestSplitZeroMismatchCell(unittest.TestCase):
    def test_empty_and_nan(self):
        self.assertEqual(_split_zero_mismatch_cell(""), [])
        self.assertEqual(_split_zero_mismatch_cell(None), [])
        self.assertEqual(_split_zero_mismatch_cell(float("nan")), [])

    def test_single_entry(self):
        cell = "TTTAGTGAAG_chrI:500:+"
        result = _split_zero_mismatch_cell(cell)
        self.assertEqual(result, ["TTTAGTGAAG"])

    def test_multiple_entries(self):
        cell = "TTTAGTGAAG_chrI:500:+,TTTGATGATG_chrX:293146:+"
        result = _split_zero_mismatch_cell(cell)
        self.assertEqual(result, ["TTTAGTGAAG", "TTTGATGATG"])

    def test_whitespace_handling(self):
        cell = " TTTAGTGAAG_chrI:500:+ , TTTGATGATG_chrX:293146:+ "
        result = _split_zero_mismatch_cell(cell)
        self.assertEqual(result, ["TTTAGTGAAG", "TTTGATGATG"])


class TestHasTtTV(unittest.TestCase):
    def test_tttv_variants(self):
        self.assertTrue(_has_ttTV(["TTTAGCG"]))  # TTTA
        self.assertTrue(_has_ttTV(["TTTCACG"]))  # TTTC
        self.assertTrue(_has_ttTV(["TTTGACG"]))  # TTTG

    def test_tttt_not_tttv(self):
        self.assertFalse(_has_ttTV(["TTTTACG"]))

    def test_empty(self):
        self.assertFalse(_has_ttTV([]))

    def test_short_seq(self):
        self.assertFalse(_has_ttTV(["TTT"]))


class TestHasTttt(unittest.TestCase):
    def test_tttt(self):
        self.assertTrue(_has_tttt(["TTTTACG"]))

    def test_tttv_not_tttt(self):
        self.assertFalse(_has_tttt(["TTTAGCG"]))

    def test_empty(self):
        self.assertFalse(_has_tttt([]))


class TestExpectedCanonicalCoordFromContig(unittest.TestCase):
    def test_plus_strand(self):
        # + strand: pos - 1
        result = _expected_canonical_coord_from_contig("chrI:501:+")
        self.assertEqual(result, "chrI:500:+")

    def test_minus_strand(self):
        # - strand: pos - 21
        result = _expected_canonical_coord_from_contig("chrIII:147715:-")
        self.assertEqual(result, "chrIII:147694:-")

    def test_invalid_contig(self):
        self.assertIsNone(_expected_canonical_coord_from_contig("garbage"))
        self.assertIsNone(_expected_canonical_coord_from_contig(""))
        self.assertIsNone(_expected_canonical_coord_from_contig(None))

    def test_pos_near_zero_plus(self):
        # pos=1, strand=+ => 1-1=0 < 1 => None
        self.assertIsNone(_expected_canonical_coord_from_contig("chrI:1:+"))

    def test_pos_near_zero_minus(self):
        # pos=10, strand=- => 10-21=-11 < 1 => None
        self.assertIsNone(_expected_canonical_coord_from_contig("chrI:10:-"))


class TestRemoveCanonicalZeroMismatch(unittest.TestCase):
    def test_removes_self_match(self):
        zcell = "TTTAGTGAAG_chrI:57761:-,TTTAGTGAAG_chrI:500:+"
        contig = "chrI:501:+"
        result = _remove_canonical_zero_mismatch(zcell, contig)
        self.assertEqual(result, "TTTAGTGAAG_chrI:57761:-")

    def test_keeps_all_when_no_match(self):
        zcell = "TTTAGTGAAG_chrI:57761:-"
        contig = "chrI:501:+"
        result = _remove_canonical_zero_mismatch(zcell, contig)
        self.assertEqual(result, "TTTAGTGAAG_chrI:57761:-")

    def test_empty_input(self):
        self.assertEqual(_remove_canonical_zero_mismatch("", "chrI:501:+"), "")
        self.assertEqual(_remove_canonical_zero_mismatch("", ""), "")


class TestCountMismatches(unittest.TestCase):
    def test_identical(self):
        seq = "TTTAGTGAAGCTTCTAGATATTTGGCG"  # 27-mer
        self.assertEqual(count_mismatches(seq, seq), 0)

    def test_known_mismatches(self):
        # After trimming 4 bases, compare up to 20 positions
        # a trim4: "AAAAAAAAAAAAAAAAAAAAAA" (22 chars), compare first 20
        # b trim4: "AAAAAAAAAAAAAAAAAAGGG" (20 chars from 24-4=20)
        a = "TTTTAAAAAAAAAAAAAAAAAAAAAA"  # 26 bases => 22 after trim
        b = "TTTTAAAAAAAAAAAAAAAAAAGAAA"  # pos 18 (0-indexed) after trim differs
        # trim4 => compare 20 chars: position 17 differs (A vs G)
        self.assertEqual(count_mismatches(a, b), 1)

    def test_multiple_mismatches(self):
        # Mismatches within the first 20 positions after trim
        a = "TTTTACGTACGTACGTACGTACGTAC"  # 25 chars
        b = "TTTTGCGTACGTACGTACGTACGTAC"  # pos 0 after trim: A→G
        self.assertEqual(count_mismatches(a, b), 1)

    def test_empty_returns_999(self):
        self.assertEqual(count_mismatches("", "ACGT"), 999)
        self.assertEqual(count_mismatches("ACGT", ""), 999)

    def test_cap_len(self):
        # Only compare first 20 bases after trim5
        a = "XXXX" + "A" * 20 + "GGG"
        b = "XXXX" + "A" * 20 + "CCC"
        self.assertEqual(count_mismatches(a, b), 0)


class TestBucketOfftargetsByMismatch(unittest.TestCase):
    def test_empty_inputs(self):
        result = bucket_offtargets_by_mismatch("", "", "")
        self.assertEqual(result, {0: [], 1: [], 2: [], 3: [], 4: []})

    def test_known_bucketing(self):
        target = "TTTTACGTACGTACGTACGTACGTAC"  # 25 bases
        same = "TTTTACGTACGTACGTACGTACGTAC"  # 0 mismatches
        # Mismatch at position 0 after trim (A→G): within cap_len=20
        one_mm = "TTTTGCGTACGTACGTACGTACGTAC"

        loci_seqs = f"{same},{one_mm}"
        loci_coords = "chrI:100:+,chrI:200:+"

        result = bucket_offtargets_by_mismatch(target, loci_seqs, loci_coords)
        self.assertEqual(len(result[0]), 1)
        self.assertIn("chrI:100:+", result[0][0])
        self.assertEqual(len(result[1]), 1)
        self.assertIn("chrI:200:+", result[1][0])

    def test_beyond_max_bucket_ignored(self):
        target = "TTTTACGTACGTACGTACGTACGTAC"
        # 5+ mismatches
        far_off = "TTTTGGGGGGGGGGGGGGGGGGTTTTT"
        result = bucket_offtargets_by_mismatch(target, far_off, "chrI:100:+")
        total = sum(len(v) for v in result.values())
        self.assertEqual(total, 0)


class TestExtractLociFromCell(unittest.TestCase):
    def test_single_forward(self):
        cell = "ACGT_1_4<chrI:500^F>"
        result = extract_loci_from_cell(cell)
        self.assertEqual(result, "chrI:500:+")

    def test_single_reverse(self):
        cell = "ACGT_1_4<chrI:500^R>"
        result = extract_loci_from_cell(cell)
        self.assertEqual(result, "chrI:500:-")

    def test_multi_site_pipe(self):
        cell = "ACGT_2_0<chrI:100^F|chrI:200^R>"
        result = extract_loci_from_cell(cell)
        self.assertEqual(result, "chrI:100:+,chrI:200:-")

    def test_multiple_tokens(self):
        cell = "ACGT_1_4<chrI:100^F>,TGCA_1_3<chrII:200^R>"
        result = extract_loci_from_cell(cell)
        self.assertEqual(result, "chrI:100:+,chrII:200:-")

    def test_empty(self):
        self.assertEqual(extract_loci_from_cell(""), "")
        self.assertEqual(extract_loci_from_cell(None), "")

    def test_real_example(self):
        # From actual flashfry data: token with multiple pipe-separated sites
        cell = "TTTCTTGAAGATTCCAGATGTTTG_2_4<chrIII:6555340^F|chrIII:6543296^F>"
        result = extract_loci_from_cell(cell)
        self.assertEqual(result, "chrIII:6555340:+,chrIII:6543296:+")


class TestParseOfftargets(unittest.TestCase):
    def test_basic(self):
        cell = "ACGTACGT_1_4<chrI:100^F>,TGCATGCA_1_3<chrII:200^R>"
        result = parse_offtargets(cell)
        self.assertEqual(result, ["ACGTACGT", "TGCATGCA"])

    def test_empty(self):
        self.assertEqual(parse_offtargets(""), [])
        self.assertEqual(parse_offtargets(None), [])


class TestBuildScoreLookup(unittest.TestCase):
    def test_basic(self):
        df = pd.DataFrame(
            {
                "RDA": ["enCas12a", "enCas12a"],
                "Pos": [1, 2],
                "MM": ["rA:dC", "rG:dT"],
                "avg_percent_active": [0.5, 0.3],
            }
        )
        name, lookup = build_score_lookup(df)
        self.assertEqual(name, "enCas12a")
        self.assertAlmostEqual(lookup[(1, "rA:dC")], 0.5)
        self.assertAlmostEqual(lookup[(2, "rG:dT")], 0.3)

    def test_missing_columns_raises(self):
        df = pd.DataFrame({"RDA": ["x"], "Pos": [1]})
        with self.assertRaises(ValueError):
            build_score_lookup(df)


class TestDnaToRnaBase(unittest.TestCase):
    def test_t_to_u(self):
        self.assertEqual(dna_to_rna_base("T"), "U")
        self.assertEqual(dna_to_rna_base("t"), "U")

    def test_others_unchanged(self):
        self.assertEqual(dna_to_rna_base("A"), "A")
        self.assertEqual(dna_to_rna_base("C"), "C")
        self.assertEqual(dna_to_rna_base("G"), "G")


class TestCfdScorePair(unittest.TestCase):
    def test_identical_sequences(self):
        seq = "TTTTACGTACGTACGTACGTACGTAC"
        lookup = {}
        self.assertAlmostEqual(cfd_score_pair(seq, seq, lookup), 1.0)

    def test_single_mismatch(self):
        target = "TTTTACGTACGTACGTACGTACGTAC"
        offtarg = "TTTTCCGTACGTACGTACGTACGTAC"
        # Position 1 after trim4: target A vs offtarget C => rA:dC
        lookup = {(1, "rA:dC"): 0.5}
        result = cfd_score_pair(target, offtarg, lookup)
        self.assertAlmostEqual(result, 0.5)

    def test_two_mismatches_multiply(self):
        target = "TTTTACGTACGTACGTACGTACGTAC"
        offtarg = "TTTTCCCTACGTACGTACGTACGTAC"
        # After trim4: pos 1: A→C = rA:dC, pos 3: G→C = rG:dC
        lookup = {(1, "rA:dC"): 0.5, (3, "rG:dC"): 0.4}
        result = cfd_score_pair(target, offtarg, lookup)
        self.assertAlmostEqual(result, 0.2)

    def test_missing_lookup_defaults_to_1(self):
        target = "TTTTACGTACGTACGTACGTACGTAC"
        offtarg = "TTTTCCGTACGTACGTACGTACGTAC"
        lookup = {}  # no entry for this mismatch
        result = cfd_score_pair(target, offtarg, lookup)
        self.assertAlmostEqual(result, 1.0)

    def test_short_seqs(self):
        self.assertAlmostEqual(cfd_score_pair("ACG", "ACG", {}), 1.0)


class TestFormatAndParseScores(unittest.TestCase):
    def test_round_trip(self):
        scores = [0.5, 1.0, 0.123]
        formatted = format_scores(scores)
        parsed = parse_score_list(formatted)
        self.assertEqual(len(parsed), 3)
        for a, b in zip(scores, parsed):
            self.assertAlmostEqual(a, b)

    def test_parse_empty(self):
        self.assertEqual(parse_score_list(""), [])
        self.assertEqual(parse_score_list(None), [])

    def test_parse_with_null(self):
        result = parse_score_list("0.5,Null,0.3")
        self.assertEqual(result, [0.5, 0.3])


class TestRoundingFunctions(unittest.TestCase):
    def test_round_float_str(self):
        self.assertEqual(_round_float_str("0.123456789"), "0.123457")
        self.assertEqual(_round_float_str("Null"), "Null")
        self.assertEqual(_round_float_str(""), "")
        self.assertEqual(_round_float_str(None), "")

    def test_round_csv_list(self):
        result = _round_csv_list("0.123456789,0.987654321,Null")
        self.assertEqual(result, "0.123457,0.987654,Null")

    def test_round_csv_empty(self):
        self.assertEqual(_round_csv_list(""), "")
        self.assertEqual(_round_csv_list(None), "")


# ---------------------------------------------------------------------------
# Genome-dependent tests (synthetic FASTA)
# ---------------------------------------------------------------------------


class TestLoadGenomeBiopython(unittest.TestCase):
    def test_loads_fasta(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
            f.write(">chrTest\nACGTACGTACGTACGTACGTACGTACGTACGT\n")
            f.flush()
            genome = load_genome_biopython(f.name)
        os.unlink(f.name)
        self.assertIn("chrTest", genome)
        self.assertEqual(len(genome["chrTest"].seq), 32)


class TestFetch27mer(unittest.TestCase):
    def setUp(self):
        # 40-base sequence for easy testing (1-based positions 1-40)
        self.genome = _make_genome_dict(
            {
                "chrT": "AAAA" + "ACGTACGTACGTACGTACGTACGTACG" + "TTTTTTTTT"
                #        1234   5                            31  32-40
            }
        )

    def test_plus_strand(self):
        # pos_1based=4, strand='+': start=5, end=31 => 27 bases
        mer = fetch_27mer(self.genome, "chrT", 4, "+")
        self.assertIsNotNone(mer)
        self.assertEqual(len(mer), 27)
        self.assertEqual(mer, "ACGTACGTACGTACGTACGTACGTACG")

    def test_minus_strand(self):
        mer = fetch_27mer(self.genome, "chrT", 4, "-")
        self.assertIsNotNone(mer)
        self.assertEqual(len(mer), 27)

    def test_missing_chrom(self):
        self.assertIsNone(fetch_27mer(self.genome, "chrZ", 4, "+"))

    def test_out_of_bounds(self):
        # pos=100 is well beyond our 40-char sequence
        self.assertIsNone(fetch_27mer(self.genome, "chrT", 100, "+"))
        # minus strand: start = pos-2 = -2 < 1 => None
        self.assertIsNone(fetch_27mer(self.genome, "chrT", 1, "-"))


class TestExpandTargetFromContig(unittest.TestCase):
    def setUp(self):
        # 60-base sequence
        self.genome = _make_genome_dict({"chrT": "N" * 5 + "ACGTACGTACGTACGTACGTACGTACG" + "N" * 28})

    def test_plus_strand(self):
        # contig 'chrT:6:+' => pos=6, strand='+' => pos-=1 => pos=5 => fetch_27mer(chrT,5,+)
        # start=6, end=32 => seq[5:32] = "ACGTACGTACGTACGTACGTACGTACG"
        result = expand_target_from_contig(self.genome, "chrT:6:+")
        self.assertEqual(len(result), 27)
        self.assertEqual(result, "ACGTACGTACGTACGTACGTACGTACG")

    def test_invalid_contig(self):
        self.assertEqual(expand_target_from_contig(self.genome, "garbage"), "")
        self.assertEqual(expand_target_from_contig(self.genome, ""), "")
        self.assertEqual(expand_target_from_contig(self.genome, None), "")


# ---------------------------------------------------------------------------
# Tests using real chrIII FASTA
# ---------------------------------------------------------------------------


@unittest.skipUnless(os.path.exists(CHRIII_FASTA), "chrIII FASTA not found")
class TestFetch27merChrIII(unittest.TestCase):
    """Validate fetch_27mer against the real ce11 chrIII sequence."""

    @classmethod
    def setUpClass(cls):
        cls.genome = load_genome_biopython(CHRIII_FASTA)

    def test_chrIII_loaded(self):
        self.assertIn("chrIII", self.genome)
        self.assertGreater(len(self.genome["chrIII"].seq), 13_000_000)

    def test_forward_fetch_known_position(self):
        # The very start of chrIII is "cctaagcctaag..."
        mer = fetch_27mer(self.genome, "chrIII", 1, "+")
        self.assertIsNotNone(mer)
        self.assertEqual(len(mer), 27)
        # positions 2-28 (0-based 1-27)
        expected = str(self.genome["chrIII"].seq[1:28]).upper()
        self.assertEqual(mer, expected)

    def test_reverse_fetch(self):
        mer = fetch_27mer(self.genome, "chrIII", 100, "-")
        self.assertIsNotNone(mer)
        self.assertEqual(len(mer), 27)


@unittest.skipUnless(os.path.exists(CHRIII_FASTA), "chrIII FASTA not found")
class TestExpandTargetChrIII(unittest.TestCase):
    """Test expand_target_from_contig using real chrIII data."""

    @classmethod
    def setUpClass(cls):
        cls.genome = load_genome_biopython(CHRIII_FASTA)

    def test_expand_known_contig(self):
        # chrIII:147715:- is a real contig from the test data
        result = expand_target_from_contig(self.genome, "chrIII:147715:-")
        self.assertEqual(len(result), 27)
        # Should be uppercase DNA
        self.assertTrue(all(c in "ACGT" for c in result))


# ---------------------------------------------------------------------------
# Integration: main() pipeline with real chrIII data
# ---------------------------------------------------------------------------


@unittest.skipUnless(
    os.path.exists(CHRIII_FASTA) and os.path.exists(MATRIX_CSV),
    "chrIII FASTA or scoring matrix not found",
)
class TestMainPipelineSynthetic(unittest.TestCase):
    """
    Run the scoring pipeline on a small synthetic FlashFry TSV that
    references only chrIII loci, using the real chrIII genome and matrix.
    """

    @classmethod
    def setUpClass(cls):
        cls.genome = load_genome_biopython(CHRIII_FASTA)
        cls.score_df = pd.read_csv(MATRIX_CSV)
        cls.matrix_name, cls.lookup = build_score_lookup(cls.score_df)

    def _make_input_tsv(self, rows: list[dict]) -> str:
        """Write a FlashFry-like TSV to a temp file and return path."""
        cols = ["contig", "start", "stop", "target", "context", "overflow", "orientation", "otCount", "offTargets"]
        df = pd.DataFrame(rows, columns=cols)
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
        df.to_csv(f, sep="\t", index=False)
        f.close()
        return f.name

    def test_rvs_rows_filtered(self):
        """RVS orientation rows should be excluded from output."""
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        rows = [
            {
                "contig": "chrIII:147715:-",
                "start": "0",
                "stop": "24",
                "target": "TTTTATGAAGCTTCAATATATTTT",
                "context": "NONE",
                "overflow": "OK",
                "orientation": "RVS",
                "otCount": "0",
                "offTargets": "",
            },
        ]
        inp = self._make_input_tsv(rows)
        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = ["prog", "-i", inp, "-m", MATRIX_CSV, "-o", out, "-g", CHRIII_FASTA]
            main()
            sys.argv = old_argv
            df = pd.read_csv(out, sep="\t")
            self.assertEqual(len(df), 0)
        finally:
            os.unlink(inp)
            os.unlink(out)

    def test_non_standard_start_stop_filtered(self):
        """Rows with start!=0 or stop!=24 should be excluded."""
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        rows = [
            {
                "contig": "chrIII:147715:-",
                "start": "1",
                "stop": "25",
                "target": "TTTTATGAAGCTTCAATATATTTT",
                "context": "NONE",
                "overflow": "OK",
                "orientation": "FWD",
                "otCount": "0",
                "offTargets": "",
            },
        ]
        inp = self._make_input_tsv(rows)
        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = ["prog", "-i", inp, "-m", MATRIX_CSV, "-o", out, "-g", CHRIII_FASTA]
            main()
            sys.argv = old_argv
            df = pd.read_csv(out, sep="\t")
            self.assertEqual(len(df), 0)
        finally:
            os.unlink(inp)
            os.unlink(out)

    def test_output_columns(self):
        """Output should have expected columns and drop offTargets."""
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        rows = [
            {
                "contig": "chrIII:147715:-",
                "start": "0",
                "stop": "24",
                "target": "TTTTATGAAGCTTCAATATATTTT",
                "context": "NONE",
                "overflow": "OK",
                "orientation": "FWD",
                "otCount": "1",
                "offTargets": "TTTTATGAAGCTTCAATATATTTT_1_0<chrIII:147715^R>",
            },
        ]
        inp = self._make_input_tsv(rows)
        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = ["prog", "-i", inp, "-m", MATRIX_CSV, "-o", out, "-g", CHRIII_FASTA]
            main()
            sys.argv = old_argv
            df = pd.read_csv(out, sep="\t")

            expected_cols = {
                "contig",
                "target",
                "otCount",
                "offTargets_loci",
                "offTargets_loci_seq",
                "0-mismatch",
                "1-mismatch",
                "2-mismatch",
                "3-mismatch",
                "4-mismatch",
                "TTTN_enCas12a",
                "TTTV_enCas12a",
                "TTTN_enCas12a_aggregated_score",
                "TTTV_enCas12a_aggregated_score",
                "unique-TTTV",
                "unique-TTTN",
            }
            self.assertTrue(
                expected_cols.issubset(set(df.columns)), f"Missing columns: {expected_cols - set(df.columns)}"
            )
            self.assertNotIn("offTargets", df.columns)
            self.assertNotIn("start", df.columns)
            self.assertNotIn("stop", df.columns)
            self.assertNotIn("orientation", df.columns)
        finally:
            os.unlink(inp)
            os.unlink(out)

    def test_canonical_self_removed_from_zero_mismatch(self):
        """The on-target locus itself should be removed from 0-mismatch."""
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        # A guide at chrIII:147715:- with itself as the only 0-mismatch hit
        rows = [
            {
                "contig": "chrIII:147715:-",
                "start": "0",
                "stop": "24",
                "target": "TTTTATGAAGCTTCAATATATTTT",
                "context": "NONE",
                "overflow": "OK",
                "orientation": "FWD",
                "otCount": "1",
                "offTargets": "TTTTATGAAGCTTCAATATATTTT_1_0<chrIII:147715^R>",
            },
        ]
        inp = self._make_input_tsv(rows)
        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = ["prog", "-i", inp, "-m", MATRIX_CSV, "-o", out, "-g", CHRIII_FASTA]
            main()
            sys.argv = old_argv
            df = pd.read_csv(out, sep="\t", dtype=str)
            zero_mm = df["0-mismatch"].iloc[0]
            # After removing canonical self-match, should be empty
            self.assertTrue(
                pd.isna(zero_mm) or str(zero_mm).strip() == "",
                f"Expected empty 0-mismatch after self-removal, got: {zero_mm}",
            )
        finally:
            os.unlink(inp)
            os.unlink(out)

    def test_unique_when_no_zero_mismatch(self):
        """A guide with no 0-mismatch hits (after self-removal) should be unique."""
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        rows = [
            {
                "contig": "chrIII:147715:-",
                "start": "0",
                "stop": "24",
                "target": "TTTTATGAAGCTTCAATATATTTT",
                "context": "NONE",
                "overflow": "OK",
                "orientation": "FWD",
                "otCount": "1",
                # Only off-target is itself (will be removed)
                "offTargets": "TTTTATGAAGCTTCAATATATTTT_1_0<chrIII:147715^R>",
            },
        ]
        inp = self._make_input_tsv(rows)
        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = ["prog", "-i", inp, "-m", MATRIX_CSV, "-o", out, "-g", CHRIII_FASTA]
            main()
            sys.argv = old_argv
            df = pd.read_csv(out, sep="\t")
            self.assertTrue(df["unique-TTTV"].iloc[0])
            self.assertTrue(df["unique-TTTN"].iloc[0])
        finally:
            os.unlink(inp)
            os.unlink(out)

    def test_aggregated_score_is_reciprocal_sum(self):
        """Aggregated score should equal 1/sum(individual scores)."""
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        # Use two chrIII off-target loci so we can verify the math
        rows = [
            {
                "contig": "chrIII:147715:-",
                "start": "0",
                "stop": "24",
                "target": "TTTTATGAAGCTTCAATATATTTT",
                "context": "NONE",
                "overflow": "OK",
                "orientation": "FWD",
                "otCount": "2",
                "offTargets": (
                    "TTTTATGAAGCTTCAATATATTTT_1_0<chrIII:147715^R>,TTTAATAATTTTTTGAATATTGGAAAA_1_4<chrIII:13211459^R>"
                ),
            },
        ]
        inp = self._make_input_tsv(rows)
        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = ["prog", "-i", inp, "-m", MATRIX_CSV, "-o", out, "-g", CHRIII_FASTA]
            main()
            sys.argv = old_argv
            df = pd.read_csv(out, sep="\t", dtype=str)

            scores_str = df["TTTN_enCas12a"].iloc[0]
            agg_str = df["TTTN_enCas12a_aggregated_score"].iloc[0]

            if scores_str and str(scores_str).strip():
                individual = parse_score_list(scores_str)
                if individual:
                    expected_agg = 1.0 / sum(individual)
                    actual_agg = float(agg_str)
                    # Aggregated score is computed from the raw (pre-rounded) values,
                    # then itself rounded to 6 decimals, so allow small tolerance
                    self.assertAlmostEqual(actual_agg, expected_agg, places=3)
        finally:
            os.unlink(inp)
            os.unlink(out)


class TestUniquenessLogic(unittest.TestCase):
    """Test unique-TTTV / unique-TTTN determination from 0-mismatch cell contents."""

    def test_no_zero_mismatch_both_unique(self):
        seqs = _split_zero_mismatch_cell("")
        self.assertEqual(len(seqs), 0)
        # Pipeline: len(seqs) < 1 => both True

    def test_tttv_match_both_not_unique(self):
        # A TTTV match (e.g. TTTA...) makes both non-unique
        seqs = ["TTTAGCGATCGATCGATCGATCGATCG"]
        self.assertTrue(_has_ttTV(seqs))
        self.assertFalse(_has_tttt(seqs))

    def test_tttt_only_tttv_unique_tttn_not(self):
        # Only TTTT matches: TTTV is still unique, TTTN is not
        seqs = ["TTTTGCGATCGATCGATCGATCGATCG"]
        self.assertFalse(_has_ttTV(seqs))
        self.assertTrue(_has_tttt(seqs))

    def test_mixed_tttv_and_tttt(self):
        seqs = ["TTTAGCG", "TTTTGCG"]
        self.assertTrue(_has_ttTV(seqs))
        self.assertTrue(_has_tttt(seqs))
        # has_ttTV => both not unique


class TestCfdScorePairWithRealMatrix(unittest.TestCase):
    """Test CFD scoring with the actual enCas12a matrix."""

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(MATRIX_CSV):
            raise unittest.SkipTest("Scoring matrix not found")
        df = pd.read_csv(MATRIX_CSV)
        cls.matrix_name, cls.lookup = build_score_lookup(df)

    def test_identical_score_is_one(self):
        seq = "TTTTACGTACGTACGTACGTACGTACG"
        self.assertAlmostEqual(cfd_score_pair(seq, seq, self.lookup), 1.0)

    def test_score_between_zero_and_one(self):
        target = "TTTTACGTACGTACGTACGTACGTACG"
        offtarg = "TTTTCCGTACGTACGTACGTACGTACG"
        score = cfd_score_pair(target, offtarg, self.lookup)
        self.assertGreater(score, 0.0)
        self.assertLess(score, 1.0)

    def test_more_mismatches_lower_or_equal_score(self):
        target = "TTTTACGTACGTACGTACGTACGTACG"
        one_mm = "TTTTCCGTACGTACGTACGTACGTACG"  # 1 mismatch at pos 1
        two_mm = "TTTTCCGTACGTACGTACGTACGTTCG"  # 2 mismatches at pos 1 and 20
        score_1 = cfd_score_pair(target, one_mm, self.lookup)
        score_2 = cfd_score_pair(target, two_mm, self.lookup)
        # More mismatches should yield a lower or equal score
        self.assertGreaterEqual(score_1, score_2)
        # And both should be < 1.0 (at least one mismatch penalty applied)
        self.assertLess(score_1, 1.0)
        self.assertLess(score_2, 1.0)


# ---------------------------------------------------------------------------
# Integration: uniqueness determination through main() with real chrIII data
# ---------------------------------------------------------------------------


@unittest.skipUnless(
    os.path.exists(CHRIII_FASTA) and os.path.exists(MATRIX_CSV) and os.path.exists(UNIQUENESS_INPUT),
    "chrIII FASTA, scoring matrix, or uniqueness test input not found",
)
class TestUniquenessThroughMain(unittest.TestCase):
    """
    Run main() on a synthetic FlashFry file with three carefully chosen
    chrIII guides that exercise all branches of the uniqueness logic:
      Row 0: TTTV 0-mismatch hit  -> unique-TTTV=False, unique-TTTN=False
      Row 1: TTTT-only 0-mismatch -> unique-TTTV=True,  unique-TTTN=False
      Row 2: self-match only       -> unique-TTTV=True,  unique-TTTN=True
    """

    def test_uniqueness_branches(self):
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = [
                "prog",
                "-i",
                UNIQUENESS_INPUT,
                "-m",
                MATRIX_CSV,
                "-o",
                out,
                "-g",
                CHRIII_FASTA,
            ]
            main()
            sys.argv = old_argv

            df = pd.read_csv(out, sep="\t")
            self.assertEqual(len(df), 3)

            # Row 0: TTTV 0-mismatch — both not unique
            self.assertFalse(df["unique-TTTV"].iloc[0])
            self.assertFalse(df["unique-TTTN"].iloc[0])

            # Row 1: TTTT-only 0-mismatch — TTTV unique, TTTN not
            self.assertTrue(df["unique-TTTV"].iloc[1])
            self.assertFalse(df["unique-TTTN"].iloc[1])

            # Row 2: no 0-mismatch after self-removal — both unique
            self.assertTrue(df["unique-TTTV"].iloc[2])
            self.assertTrue(df["unique-TTTN"].iloc[2])
        finally:
            # Clean up both the output and the slim file
            os.unlink(out)
            slim = out + ".slim.tsv.gz"
            if os.path.exists(slim):
                os.unlink(slim)

    def test_mismatch_buckets_populated(self):
        """Verify that 0-mismatch entries exist and canonical self-match is removed."""
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = [
                "prog",
                "-i",
                UNIQUENESS_INPUT,
                "-m",
                MATRIX_CSV,
                "-o",
                out,
                "-g",
                CHRIII_FASTA,
            ]
            main()
            sys.argv = old_argv

            df = pd.read_csv(out, sep="\t", dtype=str)

            # Row 0: had 2 loci, self removed => 1 remaining 0-mismatch entry
            zero_mm_0 = str(df["0-mismatch"].iloc[0]).strip()
            self.assertTrue(
                len(zero_mm_0) > 0 and zero_mm_0 != "nan", f"Expected non-empty 0-mismatch for row 0, got: {zero_mm_0}"
            )
            self.assertIn("chrIII:3066215:+", zero_mm_0)
            # Self-match at 17625 should be removed
            self.assertNotIn("17625", zero_mm_0)

            # Row 1: TTTT duplicate, self removed => 1 remaining
            zero_mm_1 = str(df["0-mismatch"].iloc[1]).strip()
            self.assertIn("chrIII:11276745:+", zero_mm_1)

            # Row 2: only self-match, removed => empty
            zero_mm_2 = df["0-mismatch"].iloc[2]
            self.assertTrue(
                pd.isna(zero_mm_2) or str(zero_mm_2).strip() == "",
                f"Expected empty 0-mismatch for row 2, got: {zero_mm_2}",
            )
        finally:
            os.unlink(out)
            slim = out + ".slim.tsv.gz"
            if os.path.exists(slim):
                os.unlink(slim)

    def test_qi_weighting_in_scores(self):
        """Verify that scores include q_i weighting (CFD * count)."""
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = [
                "prog",
                "-i",
                UNIQUENESS_INPUT,
                "-m",
                MATRIX_CSV,
                "-o",
                out,
                "-g",
                CHRIII_FASTA,
            ]
            main()
            sys.argv = old_argv

            df = pd.read_csv(out, sep="\t", dtype=str)

            # Row 0: single token with count=2, 0 mismatches => CFD=1.0
            # Weighted score should be 1.0 * 2 = 2.0
            tttn_scores = df["TTTN_enCas12a"].iloc[0]
            score_val = float(tttn_scores)
            self.assertAlmostEqual(score_val, 2.0, msg=f"Expected 2.0 (CFD=1.0 * q_i=2), got {score_val}")

            # Aggregated: 1 / 2.0 = 0.5
            agg = float(df["TTTN_enCas12a_aggregated_score"].iloc[0])
            self.assertAlmostEqual(agg, 0.5, places=4)
        finally:
            os.unlink(out)
            slim = out + ".slim.tsv.gz"
            if os.path.exists(slim):
                os.unlink(slim)


# ---------------------------------------------------------------------------
# End-to-end: compare against existing .reduced output
# ---------------------------------------------------------------------------


@unittest.skipUnless(
    os.path.exists(CHRIII_FASTA)
    and os.path.exists(MATRIX_CSV)
    and os.path.exists(FLASHFRY_INPUT)
    and os.path.exists(FLASHFRY_REDUCED),
    "Test data files not found",
)
class TestEndToEndChrIIIOnly(unittest.TestCase):
    """
    Run main() on a subset of the real FlashFry input restricted to rows
    whose off-target loci are all on chrIII, then compare with the
    existing .reduced file.

    NOTE: Since only chrIII is available, rows referencing other chromosomes
    will produce different 27-mer lookups (empty). This test validates the
    pipeline runs without error and produces the right shape/columns.
    """

    def test_pipeline_runs_on_full_input(self):
        """Pipeline completes on the full input file (loci on missing chroms
        will produce empty sequences, which is expected)."""
        from parasol_scripts.score_flashfry_cfd import main
        import sys

        out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name
        try:
            old_argv = sys.argv
            sys.argv = [
                "prog",
                "-i",
                FLASHFRY_INPUT,
                "-m",
                MATRIX_CSV,
                "-o",
                out,
                "-g",
                CHRIII_FASTA,
            ]
            main()
            sys.argv = old_argv

            df_out = pd.read_csv(out, sep="\t", dtype=str)
            df_ref = pd.read_csv(FLASHFRY_REDUCED, sep="\t", dtype=str)

            # Same number of rows (67 FWD start=0 stop=24 rows)
            self.assertEqual(len(df_out), len(df_ref))

            # Same columns
            self.assertEqual(set(df_out.columns), set(df_ref.columns))

            # Contig column should match exactly
            self.assertTrue(
                (df_out["contig"] == df_ref["contig"]).all(), "Contig columns differ between output and reference"
            )

            # otCount should match
            self.assertTrue((df_out["otCount"] == df_ref["otCount"]).all(), "otCount columns differ")

            # offTargets_loci should match (parsing is genome-independent)
            self.assertTrue(
                (df_out["offTargets_loci"].fillna("") == df_ref["offTargets_loci"].fillna("")).all(),
                "offTargets_loci columns differ",
            )
        finally:
            os.unlink(out)


if __name__ == "__main__":
    unittest.main()
