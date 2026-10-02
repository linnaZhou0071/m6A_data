"""Small fixtures for parsing, threshold/plot boundaries and MANE ID handling."""
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from filter_non_mane_nonsynonymous import clean_transcript_id
from plot_m6a_distribution import histogram
from prepare_pegg_input import parse_pos
from rebuild_source_unions import build_rows


class PipelineInvariantTests(unittest.TestCase):
    def test_plus_and_minus_site_parsing(self):
        self.assertEqual(parse_pos("chr10_100151706_-"), ("10", 100151706, "-"))
        self.assertEqual(parse_pos("chrX_42_+"), ("X", 42, "+"))

    def test_chrM_is_explicitly_unsupported_by_converter(self):
        self.assertEqual(parse_pos("chrM_16559_+"), (None, None, None))
        self.assertEqual(parse_pos("chr1_0"), (None, None, None))

    def test_mane_transcript_copy_suffix(self):
        self.assertEqual(clean_transcript_id("NM_001305203.2_4"), "NM_001305203.2")
        self.assertEqual(clean_transcript_id("NM_015299.3"), "NM_015299.3")

    def test_histogram_boundaries_include_100(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "small.csv"
            pd.DataFrame({
                "pos": ["chr1_1_+", "chr1_2_+", "chr1_3_+", "chr1_4_+"],
                "m6a_mean": [0, 9.999, 10, 100],
            }).to_csv(path, index=False)
            counts, total = histogram(path)
            self.assertEqual(total, 4)
            self.assertEqual(counts.tolist(), [2, 1, 0, 0, 0, 0, 0, 0, 0, 1])

    def test_exactly_80_from_two_methods_is_retained(self):
        rows = build_rows({
            "m6a_file1_s23": {"chr1_42_+": Decimal("70")},
            "m6a_file2": {"chr1_42_+": Decimal("90")},
        })
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["m6a_mean"], "80")
        self.assertEqual(rows[0]["in_file1_s23"], "TRUE")
        self.assertGreaterEqual(Decimal(rows[0]["m6a_mean"]), 80)

    def test_duplicate_site_fails_histogram(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicates.csv"
            pd.DataFrame({"pos": ["chr1_1_+", "chr1_1_+"], "m6a_mean": [80, 90]}).to_csv(path, index=False)
            with self.assertRaises(ValueError):
                histogram(path)


if __name__ == "__main__":
    unittest.main()
