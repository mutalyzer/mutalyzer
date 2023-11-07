import pytest

from mutalyzer.normalizer import get_position_type

EXONS = [(5026, 5113), (6010, 6127), (7020, 7165), (12958, 13948)]

@pytest.mark.parametrize(
    "position, exons, len_ss, position_type",
    [
        # in intron
        (7015, EXONS, 4, (4, 0)),
        (7016, EXONS, 4, (4, -4)),
        (7017, EXONS, 4, (4, -3)),
        (7018, EXONS, 4, (4, -2)),
        (7019, EXONS, 4, (4, -1)),
        # in exon
        (7020, EXONS, 2, (5, 1)),
        (7021, EXONS, 2, (5, 2)),
        (7022, EXONS, 2, (5, 0)),
        (7162, EXONS, 2, (5, 0)),
        (7163, EXONS, 2, (5, -2)),
        (7164, EXONS, 2, (5, -1)),
        # in intron
        (7165, EXONS, 4, (6, 1)),
        (7166, EXONS, 4, (6, 2)),
        (7167, EXONS, 4, (6, 3)),
        (7168, EXONS, 4, (6, 4)),
        (7169, EXONS, 4, (6, 0)),
    ],
)
def test_get_position_type(position, exons, len_ss, position_type):
    assert get_position_type(position, exons, len_ss) == position_type
