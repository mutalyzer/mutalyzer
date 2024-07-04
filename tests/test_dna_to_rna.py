import pytest

from mutalyzer.rna import dna_to_rna

from .commons import code_in, monkey_patches
from .variants_set import TESTS_ALL


def get_tests():
    output = []
    for test in TESTS_ALL:
        if test.get("to_test"):
            if test.get("rna_description_alt") is not None:
                if test.get("rna_description_alt") is False:
                    output.append((test["input"], None))
                else:
                    output.append((test["input"], test.get("rna_description_alt")))
            elif test.get("rna_description"):
                output.append((test["input"], test["rna_description"]))
    return output


@pytest.mark.parametrize(
    "input_description, rna_expected",
    get_tests(),
)
def test_rna(input_description, rna_expected):
    rna = dna_to_rna(input_description)
    assert rna.get("description") == rna_expected


@pytest.mark.parametrize(
    "input_description, rna_expected",
    [
        ("NG_012337.1(NM_012459.2):c.271del", "NG_012337.1(NM_012459.2):r.(269_271c[2])"),
        ("NM_012459.2:c.271del", "NM_012459.2:r.(269_271c[2])"),
        ("NG_012337.1(NM_012459.2):c.271C>A", "NG_012337.1(NM_012459.2):r.(271c>a)"),
        ("NG_012337.3(NM_003002.4):c.274G>T", "NG_012337.3(NM_003002.4):r.(274g>u)"),
        # ("NG_012337.3(NM_003002.4):c.53-20_169+10del", "NG_012337.3(NM_003002.4):r.(55_171del)"),
        ("NM_003002.4:c.274G>T", "NM_003002.4:r.(274g>u)"),
        ("NG_012337.1(NM_003002.2):c.166C>A", "NG_012337.1(NM_003002.2):r.(166c>a)"),
        # ("NG_008835.1(NM_022153.2):c.677-20_704+62del", "NG_008835.1(NM_022153.2):r.(677_704del)"),
        # ("NG_008835.1(NM_022153.2):c.83-20_511+62del", "NG_008835.1(NM_022153.2):r.(84_512del)"),
        # ("NG_008835.1(NM_022153.2):c.512-20_568+62del", "NG_008835.1(NM_022153.2):r.(512_568del)"),
        ("NR_038420.1:n.206_210del", "NR_038420.1:r.(206_210del)"),
        ("NG_007485.1(NR_024274.1):n.211delinsGG", "NG_007485.1(NR_024274.1):r.(211_215g[6])"),
    ],
)
def test_rna_new(input_description, rna_expected):
    rna = dna_to_rna(input_description)
    assert rna.get("description") == rna_expected


@pytest.mark.parametrize(
    "description",
    [
        "NG_012337.1(NM_003002.2):c.169del",
        "NG_012337.1(NM_003002.2):c.168del",
    ],
)
def test_rna_new_errors(description):
    assert dna_to_rna(description).get("errors") is not None
