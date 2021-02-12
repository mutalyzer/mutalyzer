import pytest

from normalizer.name_checker import name_check

from .commons import patch_retriever
from .variants_set import TESTS_ALL


def get_tests(tests):
    output = []
    for test in tests:
        if test.get("to_test") and test.get("coding_protein_descriptions"):
            output.append((test["input"], test["coding_protein_descriptions"]))
    return output


@pytest.mark.parametrize(
    "input_description, coding_protein_descriptions", get_tests(TESTS_ALL)
)
def test_normalizer(input_description, coding_protein_descriptions):

    normalized_output = name_check(input_description)
    normalizer_descriptions = set(normalized_output["equivalent_descriptions"]["c"])

    assert coding_protein_descriptions.issubset(normalizer_descriptions)
