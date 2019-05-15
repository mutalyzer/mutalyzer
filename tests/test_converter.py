import pytest
from .commons import get_variants_tuple, get_input_expected_tests, VARIANTS,\
 SEQUENCES
from normalizer.converter import convert_indexing, to_delins, to_hgvs


HGVS_TO_INTERNAL_AND_VICE_VERSA = get_input_expected_tests(
    [('4A>T', '3_4A>T'),
     ('4del', '3_4del'),
     ('4_5del', '3_5del'),
     ('4delA', '3_4delA'),
     ('13_16del4', '12_16del4'),
     ('4dup', '3_4dup'),
     ('3_4dup', '2_4dup'),
     ('4_5insT', '4_4insT'),
     ('4_5ins7_8', '4_4ins6_8'),
     ('4_5ins[7_8;10_20]', '4_4ins[6_8;9_20]'),
     ('4_5ins[7_8;10_20;T]', '4_4ins[6_8;9_20;T]'),
     ('11_13inv', '10_13inv'),
     ('4_5con7_8', '3_5con6_8'),
     ('4delins7_8', '3_4delins6_8'),
     ('4=', '3_4='),
     ('3_4=', '2_4=')])


@pytest.mark.parametrize('hgvs, internal', HGVS_TO_INTERNAL_AND_VICE_VERSA)
def test_to_internal_indexing(hgvs, internal):
    assert convert_indexing(hgvs, indexing='internal') == internal


@pytest.mark.parametrize('hgvs, internal', HGVS_TO_INTERNAL_AND_VICE_VERSA)
def test_to_hgvs_indexing(hgvs, internal):
    assert convert_indexing(internal, indexing='hgvs') == hgvs


@pytest.mark.parametrize('any_variant, delins', get_input_expected_tests(
    [('3_4A>T', '3_4delAinsT'),
     ('3_4del', '3_4delins'),
     ('3_5del', '3_5delins'),
     ('3_4dup', '4_4delins3_4'),
     ('15_15insT', '15_15delinsT'),
     ('4_4ins6_8', '4_4delins6_8'),
     ('4_4ins[6_8;9_20]', '4_4delins[6_8;9_20]'),
     ('10_13inv', '10_13delins10_13inv'),
     ('3_5con6_8', '3_5delins6_8'),
     ('3_4delins6_8', '3_4delins6_8'),
     ('3_4=', '3_4delins3_4')]))
def test_to_delins(any_variant, delins):
    assert to_delins(any_variant) == delins


@pytest.mark.parametrize('any_variant, delins', get_input_expected_tests(
    [('3_4del', '3_4delins'),
     ('3_5del', '3_5delins'),
     ('15_15insT', '15_15delinsT'),
     ('4_4ins6_8', '4_4delins6_8'),
     ('4_4ins[6_8;9_20]', '4_4delins[6_8;9_20]'),
     ('10_13inv', '10_13delins10_13inv'),
     ('3_4delins6_8', '3_4delins6_8'),
     ('3_4=', '3_4delins3_4')]))
def test_to_hgvs_no_sequences(any_variant, delins):
     assert to_hgvs(delins) == any_variant
