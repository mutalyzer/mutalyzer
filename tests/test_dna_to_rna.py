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
    [
        ("NG_012337.1(NM_012459.2):c.271del", "NG_012337.1(NM_012459.2):r.(269_271c[2])"),
        ("NM_012459.2:c.271del", "NM_012459.2:r.(269_271c[2])"),
        ("NG_012337.3(NM_003002.4):c.274G>T", "NG_012337.3(NM_003002.4):r.(274g>u)"),
        ("NM_003002.4:c.274G>T", "NM_003002.4:r.(274g>u)"),
        ("NG_012337.1(NM_003002.2):c.169T>A", "NG_012337.1(NM_003002.2):r.169u>a"),
    ],
)
def test_rna_new(input_description, rna_expected):
    rna = dna_to_rna(input_description)
    assert rna == rna_expected
