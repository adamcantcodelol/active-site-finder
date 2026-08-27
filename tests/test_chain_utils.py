from chain_utils import choose_chain, normalize_chain, validate_chain_residues


def test_choose_longest_protein_chain():
    result = choose_chain([
        {"chain": "A", "length": 40, "is_protein": True},
        {"chain": "B", "length": 120, "is_protein": True},
    ])
    assert result["chain"] == "B"


def test_preferred_chain_wins():
    result = choose_chain([
        {"chain": "A", "length": 200, "is_protein": True},
        {"chain": "B", "length": 80, "is_protein": True},
    ], "B")
    assert result["chain"] == "B"


def test_normalize_chain():
    assert normalize_chain(" A ") == "A"
    assert normalize_chain("  ") is None


def test_validate_residue_records():
    assert validate_chain_residues([{"resname": "HIS", "resnum": 57}])
    assert not validate_chain_residues([{"resname": "ATP", "resnum": 1}])
    assert not validate_chain_residues([])
