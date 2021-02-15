import pytest

from normalizer.name_checker import name_check

from .commons import code_in, patch_retriever
from .variants_set import TESTS_ALL

TESTS_ERROR = [
    ("NG_007485.1:g.0del", "EOUTOFBOUNDARY"),
    ("NG_007485.1:g.-1del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NM_000077.4):c.40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NM_000077.4):c.-40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NR_003529.3):n.40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NR_003529.3):n.-40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1:g.40_42insT", "EINSERTIONRANGE"),
    ("NG_007485.1(NM_000077.4):c.40_42insT", "EINSERTIONRANGE"),
    # Different exons positive strand.
    ("NG_007485.1(NM_000077.4):c.150_151insT", "EINSERTIONRANGE"),
    # Different exons negative strand.
    ("NG_012337.1(NM_012459.2):c.129_130insT", "EINSERTIONRANGE"),
    ("NG_007485.1(NR_003529.3):n.40_42insT", "EINSERTIONRANGE"),
    ("NG_012337.1(NM_003002.2):c.100_102AAT[4]", "EREPEATMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.100_102AA[4]", "EREPEATREFERENCELENGTH"),
    ("NG_012337.1(NM_003002.2):c.274A>T", "EDELETEDSEQUENCEMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274ATT>T", "EDELETEDSEQUENCEMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274del3", "EDELETEDLENGTHMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274_279del3", "EDELETEDLENGTHMISMATCH"),
    ("NG_012337.1:g.7125delGACinsT", "EDELETEDSEQUENCEMISMATCH"),
    ("NG_012337.1:g.[10_20del;12_13insA]", "EOVERLAP"),
    ("NG_012337.1(NM_003002.2):c.[274G>T;274del]", "EOVERLAP")
    # ("NG_007485.1:g.4_3del", "ERANGE"),
]

TESTS_NO_ERROR = [
    "NG_007485.1:g.33741del",
    "NG_007485.1:g.1del",
    "NG_007485.1(NM_000077.4):c.-24664del",
    "NG_007485.1(NM_000077.4):c.149_150insT",
    "NG_007485.1(NM_000077.4):c.151_152insT",
    "NG_012337.1(NM_012459.2):c.128_129insT",
    "NG_012337.1(NM_003002.2):c.274G>T",
    "NG_012337.1(NM_003002.2):c.274del1",
    "NG_012337.1:g.7125delGinsT",
    "NG_012337.1(NM_003002.2):c.274_275del2",
    "NG_012337.1(NM_012459.2):c.130_131insT",
]


@pytest.mark.parametrize("input_description, code", TESTS_ERROR)
def test_error(input_description, code):
    assert code_in(code, name_check(input_description)["errors"])


@pytest.mark.parametrize("input_description", TESTS_NO_ERROR)
def test_no_errors(input_description):
    assert name_check(input_description).get("errors") is None


def get_tests(tests, code_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(code_type):
            output.append((test["input"], test[code_type]))
    return output


@pytest.mark.parametrize("input_description, codes", get_tests(TESTS_ALL, "errors"))
def test_errors_from_variants_set(input_description, codes):
    for code in codes:
        assert code_in(code, name_check(input_description)["errors"])


@pytest.mark.parametrize("input_description, codes", get_tests(TESTS_ALL, "infos"))
def test_infos_from_variants_set(input_description, codes):
    for code in codes:
        assert code_in(code, name_check(input_description)["infos"])
