from chimerax_ui import molstar_url


def test_molstar_url_pdb_only():
    assert molstar_url("2QRU") == "https://molstar.org/viewer/?pdb=2qru"


def test_molstar_url_with_chain_and_residues():
    url = molstar_url("2QRU", "A", [57, 58, 59])
    assert "pdb=2qru" in url
    assert "chain=A" in url
    assert "residues=57%2C58%2C59" in url
