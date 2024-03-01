import pytest

from mutalyzer.rna import dna_to_rna

from .commons import code_in, monkey_patches
from .variants_set import TESTS_ALL


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(t_type):
            output.append((test["input"], test[t_type]))
    return output


@pytest.mark.parametrize(
    "input_description, rna_expected",
    get_tests(TESTS_ALL, "rna_description"),
)
def test_rna(input_description, rna_expected):
    rna = dna_to_rna(input_description)
    assert rna == rna_expected


@pytest.mark.parametrize(
    "input_description, rna_expected",
    [("NG_012337.1(NM_012459.2):c.271del", "NG_012337.1(NM_012459.2):r.271del")],
)
def test_rna_new(input_description, rna_expected):
    rna = dna_to_rna(input_description)
    assert rna == rna_expected
