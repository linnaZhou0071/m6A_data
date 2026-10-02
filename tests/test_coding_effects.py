"""Small strand and CDS-boundary tests for corrected A-to-G effect logic."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from predict_AtoG_coding_effects_v2 import (
    compute_coding_offset,
    genomic_bases_to_mrna_codon,
    offsets_to_genomic_coords,
)


class CodingEffectTests(unittest.TestCase):
    def test_forward_codon_is_unchanged(self):
        self.assertEqual(genomic_bases_to_mrna_codon("ACT", "+"), "ACT")

    def test_reverse_codon_is_complemented_not_reversed_twice(self):
        # Genome bases are already in transcript order at coordinates 202, 201, 200.
        # Their transcript codon is the complement, ACT; reverse-complement gives TCA.
        self.assertEqual(genomic_bases_to_mrna_codon("TGA", "-"), "ACT")

    def test_reverse_exon_boundary_coordinates(self):
        exons = [(100, 103), (200, 203)]
        self.assertEqual(compute_coding_offset(exons, 203, "-"), 0)
        self.assertEqual(compute_coding_offset(exons, 201, "-"), 2)
        self.assertEqual(compute_coding_offset(exons, 103, "-"), 3)
        self.assertEqual(offsets_to_genomic_coords(exons, "-", 1, 3), [202, 201, 103])

    def test_forward_exon_boundary_coordinates(self):
        exons = [(100, 103), (200, 203)]
        self.assertEqual(compute_coding_offset(exons, 101, "+"), 0)
        self.assertEqual(compute_coding_offset(exons, 202, "+"), 4)
        self.assertEqual(offsets_to_genomic_coords(exons, "+", 2, 3), [103, 201, 202])

    def test_invalid_strand_fails(self):
        with self.assertRaises(ValueError):
            genomic_bases_to_mrna_codon("ACT", "?")


if __name__ == "__main__":
    unittest.main()
