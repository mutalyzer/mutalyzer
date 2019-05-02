import pytest
from .commons import get_variants_tuple
from normalizer.converter import convert_indexing

HGVS_INTERNAL = [get_variants_tuple('4del', '3_4del'),
                 get_variants_tuple('4_5del', '3_5del'),
                 get_variants_tuple('4A>T', '3_4A>T'),
                 get_variants_tuple('4_5insT', '4_4insT'),
                 get_variants_tuple('4_5ins7_8', '4_4ins6_8'),
                 get_variants_tuple('4_5ins[7_8;10_20]', '4_4ins[6_8;9_20]'),
                 get_variants_tuple('4_5ins[7_8;10_20;T]',
                                    '4_4ins[6_8;9_20;T]'),
                 get_variants_tuple('4=', '3_4='),
                 get_variants_tuple('3_4=', '2_4='),
                 get_variants_tuple('4dup', '3_4dup'),
                 get_variants_tuple('3_4dup', '2_4dup'),
                 get_variants_tuple('4_5con7_8', '3_5con6_8'),
                 get_variants_tuple('4delins7_8', '3_4delins6_8')]


@pytest.mark.parametrize('hgvs, internal', HGVS_INTERNAL)
def test_to_internal_indexing(hgvs, internal):
    assert convert_indexing(hgvs, indexing='internal') == internal


@pytest.mark.parametrize('hgvs, internal', HGVS_INTERNAL)
def test_to_hgvs_indexing(hgvs, internal):
    assert convert_indexing(internal, indexing='hgvs') == hgvs
