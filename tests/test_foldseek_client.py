from foldseek_client import _alignment_start_score


def _residues(one_letter_sequence):
    return [{"one": letter} for letter in one_letter_sequence]


def test_anchor_finds_position_near_end_of_short_chain():
    # The alignment only overlaps the LAST two residues of a 3-residue chain.
    # A version of this function that requires the whole probe to fit inside
    # the residue list can never test that position and wrongly reports no
    # match at all, even though 2 of 3 letters do line up.
    residues = _residues("GHD")  # G(1) H(2) D(3)
    idx, identity = _alignment_start_score("HDA", residues, start_hint=1)
    assert idx == 1
    assert identity >= 0.60


def test_anchor_still_rejects_genuinely_unrelated_sequence():
    residues = _residues("GHD")
    idx, identity = _alignment_start_score("QQQ", residues, start_hint=1)
    assert idx == -1
    assert identity < 0.60


def test_anchor_handles_empty_inputs():
    assert _alignment_start_score("", [{"one": "A"}], 1) == (-1, 0.0)
    assert _alignment_start_score("AHD", [], 1) == (-1, 0.0)
