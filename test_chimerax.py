from chimerax import residue_spec, select_command


def test_selection_command_preserves_insertion_codes():
    assert select_command("#2", "D", [(57, "A"), (57, ""), (103, "")]) == "select #2/D:57,57A,103"


def test_selection_command_deduplicates_and_sorts():
    assert select_command("1", "A", [103, 57, 103, "57"]) == "select #1/A:57,103"


def test_empty_selection_is_safe():
    assert select_command("#1", "A", []) == "select clear"


def test_bad_model_and_chain_are_safely_normalized():
    assert select_command("not-a-model", "A/B", [12]) == "select #1:12"


def test_residue_spec_keeps_valid_insertion_code():
    assert residue_spec(42, "B") == "42B"
