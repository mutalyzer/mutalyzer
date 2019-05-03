import pytest
from .commons import get_variants_tuple
from normalizer.converter import convert_indexing, to_delins


HGVS_INTERNAL = [get_variants_tuple('4A>T', '3_4A>T'),
                 get_variants_tuple('4del', '3_4del'),
                 get_variants_tuple('4_5del', '3_5del'),
                 get_variants_tuple('4delA', '3_4delA'),
                 get_variants_tuple('13_16del4', '12_16del4'),
                 get_variants_tuple('4dup', '3_4dup'),
                 get_variants_tuple('3_4dup', '2_4dup'),
                 get_variants_tuple('4_5insT', '4_4insT'),
                 get_variants_tuple('4_5ins7_8', '4_4ins6_8'),
                 get_variants_tuple('4_5ins[7_8;10_20]', '4_4ins[6_8;9_20]'),
                 get_variants_tuple('4_5ins[7_8;10_20;T]',
                                    '4_4ins[6_8;9_20;T]'),
                 get_variants_tuple('11_13inv', '10_13inv'),
                 get_variants_tuple('4_5con7_8', '3_5con6_8'),
                 get_variants_tuple('4delins7_8', '3_4delins6_8'),
                 get_variants_tuple('4=', '3_4='),
                 get_variants_tuple('3_4=', '2_4=')]


@pytest.mark.parametrize('hgvs, internal', HGVS_INTERNAL)
def test_to_internal_indexing(hgvs, internal):
    assert convert_indexing(hgvs, indexing='internal') == internal


@pytest.mark.parametrize('hgvs, internal', HGVS_INTERNAL)
def test_to_hgvs_indexing(hgvs, internal):
    assert convert_indexing(internal, indexing='hgvs') == hgvs


TO_DELINS = [get_variants_tuple('3_4A>T', '3_4delAinsT'),
             get_variants_tuple('3_4del', '3_4delins'),
             get_variants_tuple('3_5del', '3_5delins'),
             get_variants_tuple('3_4dup', '4_4delins3_4'),
             get_variants_tuple('15_15insT', '15_15delinsT'),
             get_variants_tuple('4_5ins7_8', '4_5delins7_8'),
             get_variants_tuple('4_4ins[6_8;9_20]', '4_4delins[6_8;9_20]'),
             get_variants_tuple('10_13inv', '10_13delins10_13inv'),
             get_variants_tuple('3_5con6_8', '3_5delins6_8'),
             get_variants_tuple('3_4delins6_8', '3_4delins6_8'),
             get_variants_tuple('3_4=', '3_4delins3_4')]


@pytest.mark.parametrize('any_variant, delin', TO_DELINS)
def test_to_delins(any_variant, delin):
    assert to_delins(any_variant) == delin
