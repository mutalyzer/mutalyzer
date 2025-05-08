import pytest

from mutalyzer.normalizer import normalize_alt

from .commons import code_in, monkey_patches
from .variants_set import TESTS_ALL


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(t_type):
            output.append((test["input"], test[t_type]))
    return output


@pytest.mark.parametrize(
    "input_description, normalized_alt", get_tests(TESTS_ALL, "normalized_alt")
)
def test_normalize_alt(input_description, normalized_alt):
    d = normalize_alt(input_description)
    assert d["normalized_description"] == normalized_alt


@pytest.mark.parametrize("input_description, genomic", get_tests(TESTS_ALL, "genomic_alt"))
def test_genomic_alt(input_description, genomic):
    d = normalize_alt(input_description)
    if d["equivalent_descriptions"].get("g"):
        assert d["equivalent_descriptions"]["g"][0]["description"] == genomic
