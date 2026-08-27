from site_match import SitePair, choose_local_triplet, correspondence_is_valid, format_pair


def pair(q, h, exact=False):
    return SitePair("D", h, "", "HIS", "A", q, "", "HIS", exact)


def test_choose_local_triplet_is_exactly_three():
    pairs = [pair(10, 20), pair(11, 21), pair(12, 22), pair(100, 110)]
    chosen = choose_local_triplet(pairs)
    assert len(chosen) == 3
    assert [p.query_resnum for p in chosen] == [10, 11, 12]


def test_choose_local_triplet_prefers_local_both_structures():
    pairs = [pair(10, 20), pair(11, 21), pair(50, 22), pair(51, 23), pair(52, 24)]
    chosen = choose_local_triplet(pairs)
    assert [p.query_resnum for p in chosen] == [50, 51, 52]


def test_conservation_is_tiebreaker():
    pairs = [pair(10, 20, False), pair(11, 21, False), pair(12, 22, False), pair(13, 23, True)]
    chosen = choose_local_triplet(pairs)
    assert len(chosen) == 3
    assert chosen[-1].exact is True


def test_alignment_column_must_match_real_residues():
    assert correspondence_is_valid({"resname": "HIS"}, {"resname": "HIS"}, "H", "H")
    assert not correspondence_is_valid({"resname": "ASP"}, {"resname": "HIS"}, "H", "H")
    assert not correspondence_is_valid({"resname": "HIS"}, {"resname": "HIS"}, "-", "H")


def test_format_pair_is_sprite_style():
    text = format_pair(pair(57, 57, True))
    assert text == "D57 HIS matches A57 HIS ✓"
