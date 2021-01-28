import pytest

from normalizer.name_checker import name_check

from .commons import patch_retriever


@pytest.mark.parametrize(
    "hgvs_description, normalized_description",
    [
        ("NG_012337.1(SDHD_v001):c.274G>T", "NG_012337.1:g.7125G>T"),
        # ('NG_012337.1(NM_003002.2):c.274G>T',
        #  'NG_012337.1:g.7125G>T'),
        # Roll
        # ----
        # We should roll.
        ("NM_003002.4:c.273del", "NM_003002.4:n.309del"),
        # We cannot roll, we are at the end.
        ("NM_003002.4:c.274del", "NM_003002.4:n.309del"),
        # We can roll but should not, because it is over a splice site.
        ("NM_000088.3:c.333del", "NM_000088.3:n.459del"),
        # We can roll two positions, but should roll only one because
        # otherwise it is over a splice site.
        ("NM_000088.3:c.368del", "NM_000088.3:n.495del"),
        # We can roll and should, we stay in the same exon.
        ("NM_000088.3:c.334del", "NM_000088.3:n.461del"),
        # ('NG_012337.1:c.274G>T',
        #  'NG_012337.1:g.7125G>T'),
        # ('NC_000024.10(PRY_v001):c.948delC',
        #  'NC_000024.10:g.22514574del'),
        # ('NC_000024.10(XR_001756027.1):c.948del',
        #  'NC_000024.10:g.1352647del'),
        # ('NC_000001.11(ZRANB2-AS1_v001):c.948delC',
        #  'NC_000001.11:g.71047454del'),
        # ('NC_000001.11(NR_038420.1):c.948delC',
        #  'NC_000001.11:g.71047454del'),
        # To be added
        # ('NG_009497.1(LOC102723833_v001):c.8682-19dup',
        #  'NG_009497.1:g.561208dup'),
        # ('NM_003002.4:c.5762_5763insNG_009113.2:g.91138_91274',
        #  'NG_009497.1:g.561208dup'),
        # ('NG_009113.2:g.575_576insNG_009113.2:g.9118_9127',
        #  'NG_009497.1:g.561208dup'),
        # ('NG_009113.2:g.575_576insNG_009497.1:g.9118_9127',
        #  'NG_009497.1:g.561208dup'),
        # ('NG_017013.2(NM_001126118.1):c.100-10del',
        #  'Offset may not be from position 100 because this is not an exon boundary.'),
        # ('NG_017013.2(NM_001126118.1):c.259del',
        #  'Offset may not be from position 100 because this is not an exon boundary.'),
    ],
)
def test_mutalyzer3(hgvs_description, normalized_description):
    print(name_check(hgvs_description))
    assert name_check(hgvs_description) == normalized_description
