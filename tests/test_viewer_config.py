from viewer_config import molstar_focus_url, residue_query, viewer_selection_payload


def test_residue_query_deduplicates():
    assert residue_query([57, "57", "58A", "", None]) == "57,58A"


def test_molstar_focus_url():
    url = molstar_focus_url("2QRU", "A", [57, 58, "59A"])
    assert url.startswith("https://molstar.org/viewer/?")
    assert "pdb=2qru" in url
    assert "chain=A" in url
    assert "residues=57%2C58%2C59A" in url


def test_viewer_payload():
    assert viewer_selection_payload("A", [57, 58]) == {"chain": "A", "residues": ["57", "58"]}
