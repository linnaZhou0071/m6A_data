"""Smoke tests for the guide-only Zenodo export."""
import csv
import tempfile
import unittest
from pathlib import Path

from scripts.export_zenodo_candidates import FIELDS, export


class ExportZenodoCandidatesTests(unittest.TestCase):
    def make_source(self, root: Path, allele: str = "G") -> Path:
        source = root / "source.csv"
        row = {field: "1" for field in FIELDS}
        row.update({
            "Original_Pos_ID": "chr1_123_+", "pegRNA_id": "pegRNA_0",
            "Chromosome": "1", "Start_Position": "123",
            "Reference_Allele": "A", "Tumor_Seq_Allele2": allele,
            "Protospacer": "ACGT", "PAM": "AGG", "RTT": "ACGT",
            "PBS": "CGTA", "wt_w_context": "MUST_NOT_EXPORT",
        })
        with source.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(FIELDS) + ["wt_w_context"])
            writer.writeheader()
            writer.writerow(row)
        return source

    def test_export_is_lean_and_counted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = export(self.make_source(root), root / "out.csv", 1, 1)
            self.assertEqual((result["rows"], result["sites"]), (1, 1))
            with (root / "out.csv").open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(tuple(reader.fieldnames), FIELDS)
                self.assertEqual(next(reader)["Original_Pos_ID"], "chr1_123_+")

    def test_bad_allele_does_not_publish_partial_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "out.csv"
            with self.assertRaisesRegex(ValueError, "Allele mismatch"):
                export(self.make_source(root, "T"), target, 1, 1)
            self.assertFalse(target.exists())


if __name__ == "__main__":
    unittest.main()
