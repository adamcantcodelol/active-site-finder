from structure_viewer import normalize_pdb_id, normalize_residues, viewer_html


def test_normalize_pdb_id():
    assert normalize_pdb_id("2QRU") == "2qru"


def test_normalize_residues_deduplicates():
    assert normalize_residues([57, "58A", 57, "bad!"]) == ["57", "58A"]


def test_viewer_contains_selected_chain_and_residues():
    html = viewer_html("2QRU", "A", [57, 58, 59])
    assert "3Dmol-min.js" in html
    assert '"A"' in html
    assert '"57"' in html
    assert '"58"' in html
    assert '"59"' in html
    assert "zoomTo" in html
