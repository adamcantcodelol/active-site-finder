import unittest

from active_site import _anchor_alignment, map_binding_site_details, select_local_triplet
from foldseek_client import Hit


class ActiveSiteMappingTests(unittest.TestCase):
    def test_anchor_handles_unresolved_residues_before_alignment(self):
        residues = [
            {"resnum": 10, "insertion_code": "", "resname": "GLY", "one": "G"},
            {"resnum": 11, "insertion_code": "", "resname": "ALA", "one": "A"},
            {"resnum": 12, "insertion_code": "", "resname": "HIS", "one": "H"},
            {"resnum": 13, "insertion_code": "", "resname": "ASP", "one": "D"},
        ]
        idx, confidence = _anchor_alignment("AHD", residues, 1)
        self.assertEqual(idx, 1)
        self.assertGreater(confidence, 0.9)

    def test_mapping_rejects_shifted_sequence(self):
        target = """ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  CA  HIS A   2       1.000   0.000   0.000  1.00 20.00           C
ATOM      3  CA  ASP A   3       2.000   0.000   0.000  1.00 20.00           C
SITE     1 AC1  2 HIS A   2  ASP A   3
"""
        query = """ATOM      1  CA  ALA B  20       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  CA  HIS B  21       1.000   0.000   0.000  1.00 20.00           C
ATOM      3  CA  ASP B  22       2.000   0.000   0.000  1.00 20.00           C
"""
        hit = Hit("TEST_A", "test", "pdb100", None, 1.0, None, None, 1.0, 1, 3, 1, 3, "HDA", "HDA")
        mapped = map_binding_site_details(hit, target, query, "B")
        self.assertEqual([(p["tn"], p["qn"]) for p in mapped], [(2, 21), (3, 22)])
        self.assertTrue(all(p["exact"] for p in mapped))

    @staticmethod
    def _pair(qn, tn, exact=False):
        return {"qn": qn, "qicode": "", "qname": "ALA", "tn": tn,
                "ticode": "", "tname": "VAL", "exact": exact}

    def test_local_selection_is_exactly_three_residues(self):
        pairs = [self._pair(10, 50), self._pair(11, 51), self._pair(12, 52), self._pair(80, 120)]
        result = select_local_triplet(pairs)
        self.assertEqual(len(result), 3)
        self.assertEqual([p["qn"] for p in result], [10, 11, 12])

    def test_local_selection_prefers_closeness_over_distant_conservation(self):
        pairs = [self._pair(10, 50, True), self._pair(11, 51, True),
                 self._pair(12, 52, True), self._pair(200, 250, True)]
        result = select_local_triplet(pairs)
        self.assertEqual([p["qn"] for p in result], [10, 11, 12])

    def test_empty_local_site_is_safe(self):
        self.assertEqual(select_local_triplet([]), [])


if __name__ == "__main__":
    unittest.main()
