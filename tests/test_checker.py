import pytest
from .commons import get_variant
from normalizer.checker import is_overlap, are_sorted


def get_variants_tests(tests):
    output = []
    for test in tests:
        test_list = []
        for variant_key in test[0]:
            test_list.append(get_variant(variant_key))
        output.append((test_list, test[1]))
    return output


@pytest.mark.parametrize('variants, expected_result',
                         get_variants_tests(
                             [(['4_4delinsAA', '10del'], True),
                              (['10del', '(10_20)_(30_40)del'], True),
                              (['10del', '4_4delinsAA'], False),
                              (['10del', '4_4delinsAA'], False)]))
def test_is_sorted(variants, expected_result):
    assert are_sorted(variants) == expected_result


@pytest.mark.parametrize('variants, expected_result',
                         get_variants_tests(
                             [(['4_4delinsAA', '10del'], False),
                              (['(10_15)del', '10del'], True),
                              (['10_20=', '10del'], True)]))
def test_is_overlap(variants, expected_result):
    assert is_overlap(variants) == expected_result


def test_semantic_check():
    pass


def test_reference_check():
    pass
