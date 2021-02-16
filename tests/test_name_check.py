import pytest

from normalizer.name_checker import name_check

from .commons import code_in, patch_retriever
from .variants_set import TESTS_ALL


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(t_type):
            output.append((test["input"], test[t_type]))
    return output


@pytest.mark.parametrize(
    "input_description, normalized", get_tests(TESTS_ALL, "normalized")
)
def test_normalize(input_description, normalized):
    d = name_check(input_description)
    assert d["normalized_description"] == normalized


@pytest.mark.parametrize("input_description, genomic", get_tests(TESTS_ALL, "genomic"))
def test_genomic(input_description, genomic):
    d = name_check(input_description)
    if d["equivalent_descriptions"].get("g"):
        assert d["equivalent_descriptions"]["g"][0] == genomic


@pytest.mark.parametrize("input_description, genomic", get_tests(TESTS_ALL, "genomic"))
def test_genomic(input_description, genomic):
    d = name_check(input_description)
    if d["equivalent_descriptions"].get("g"):
        assert d["equivalent_descriptions"]["g"][0] == genomic


@pytest.mark.parametrize("input_description, codes", get_tests(TESTS_ALL, "errors"))
def test_errors(input_description, codes):
    for code in codes:
        assert code_in(code, name_check(input_description)["errors"])


@pytest.mark.parametrize("input_description, codes", get_tests(TESTS_ALL, "infos"))
def test_infos(input_description, codes):
    for code in codes:
        assert code_in(code, name_check(input_description)["infos"])
