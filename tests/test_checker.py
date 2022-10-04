import pytest

from mutalyzer.normalizer import normalize

from .commons import code_in, patch_retriever
from .variants_set import TESTS_ALL

TESTS_ERROR = [
    ("NG_007485.1:g.0del", "EOUTOFBOUNDARY"),
    ("NG_007485.1:g.-1del", "EOUTSIDECDS"),
    ("NG_007485.1:g.1000000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NM_000077.4):c.40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NM_000077.4):c.-40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NR_003529.3):n.40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NR_003529.3):n.-40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NR_003529.3):n.10delinsNG_007485.1:g.0", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NR_003529.3):n.10delinsNG_007485.1:g.0", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NR_003529.3):n.10delinsNG_007485.1:g.100000", "EOUTOFBOUNDARY"),
    ("NG_012337.1(NM_003002.2):c.*2825T>T", "EOUTOFBOUNDARY"),
    ("NG_012337.1:g.15949del", "EOUTOFBOUNDARY"),
    ("NG_012337.1(NM_003002.2):c.-5062G>T", "EOUTOFBOUNDARY"),
    ("NG_007485.1:g.4_3del", "ERANGEREVERSED"),
    ("NG_007485.1(NM_000077.4):c.135_130insT", "ERANGEREVERSED"),
    # Negative strand.
    ("NG_012337.1(NM_012459.2):c.135_130insT", "ERANGEREVERSED"),
    (
        "NG_012337.1(NM_012459.2):c.1_2ins[LRG_1:g.100000;AAA]",
        "EOUTOFBOUNDARY",
    ),
    ("NG_012337.1(NM_012459.2):c.-11025_11020inv", "EOUTOFBOUNDARY"),
    ("NG_007485.1:g.40_42insT", "EINSERTIONRANGE"),
    ("NG_007485.1(NM_000077.4):c.40_42insT", "EINSERTIONRANGE"),
    # Different exons positive strand.
    ("NG_007485.1(NM_000077.4):c.150_151insT", "EINSERTIONRANGE"),
    # Different exons negative strand.
    ("NG_012337.1(NM_012459.2):c.129_130insT", "EINSERTIONRANGE"),
    ("NG_007485.1(NR_003529.3):n.40_42insT", "EINSERTIONRANGE"),
    ("NG_012337.1(NM_003002.2):c.100_102AAT[4]", "EREPEATMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.100_102AA[4]", "EREPEATREFERENCELENGTH"),
    ("NG_012337.1(NM_003002.2):c.274A>T", "ESEQUENCEMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274ATT>T", "ESEQUENCEMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274del3", "ELENGTHMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274_279del3", "ELENGTHMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274inv3", "ELENGTHMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274_279inv3", "ELENGTHMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274_279invAA", "ESEQUENCEMISMATCH"),
    ("NG_012337.1(NM_003002.2):c.274invAA", "ESEQUENCEMISMATCH"),
    ("NG_012337.1:g.7125delGACinsT", "ESEQUENCEMISMATCH"),
    ("NG_012337.1:g.[10_20del;12_13insA]", "EOVERLAP"),
    ("NG_012337.1(NM_003002.2):c.[274G>T;274del]", "EOVERLAP"),
    ("NG_012337.1:g.7125+1G>T", "EOFFSET"),
    ("NG_012337.1:g.5+3del", "EOFFSET"),
    ("NG_012337.1:g.4_5+3del", "EOFFSET"),
    ("NG_012337.1:g.5+3_7del", "EOFFSET"),
    ("NG_012337.1:g.10_11insLRG_24:g.2+5", "EOFFSET"),
    ("NG_012337.1:g.274dup[1;A]", "EINSERTEDLENGTH"),
    ("NG_012337.1:g.20_21ins[(123);30_40]", "EINSERTEDLENGTH"),
    ("NG_012337.1:g.20_21ins[123;30_40]", "EINSERTEDLENGTH"),
    ("NG_012337.1:g.20_21ins[30_40;123]", "EINSERTEDLENGTH"),
    ("NG_012337.1:g.20_21ins[?]", "EINSERTEDLENGTH"),
    ("NG_012337.1:g.20_21ins[30_40;?]", "EINSERTEDLENGTH"),
    ("NG_012337.1:g.20>60", "EINSERTEDLENGTH"),
    ("NG_012337.1:g.20_21ins[30_?]", "EUNCERTAIN"),
    ("NG_012337.1:g.20_21ins[?_40]", "EUNCERTAIN"),
    ("NG_012337.1:g.20_21ins[30_40;50_?]", "EUNCERTAIN"),
    ("NG_012337.1:g.20_21ins[30_40;?_60]", "EUNCERTAIN"),
    ("NM_004152.3:c.205_685del", "ECDSSLICES"),
    ("NM_004152.3:p.205_685del", "ECDSSLICES"),
]

TESTS_NO_ERROR = [
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
    "NG_012337.1(NM_003002.2):c.*2824T>T",
    "NG_012337.1:g.15948del",
    "NG_012337.1(NM_003002.2):c.-5061G>T",
    "NG_012337.1:g.0_1insAAA",
    "NG_012337.1:g.15948_15949insAAA",
    "NG_012337.1:g.274dup1",
    "NG_012337.1:g.274dupT",
    "NG_012337.1:g.274inv1",
    "NG_012337.1:g.274invT",
    "NG_012337.1:g.274_276inv3",
    "NG_012337.1:g.274_276invTAC",
    "NM_004152.3:n.200del",
    "NM_004152.3:r.200del",
]


@pytest.mark.parametrize("input_description, code", TESTS_ERROR)
def test_error(input_description, code):
    assert code_in(code, normalize(input_description)["errors"])


@pytest.mark.parametrize("input_description", TESTS_NO_ERROR)
def test_no_errors(input_description):
    assert normalize(input_description).get("errors") is None


def get_tests(tests, code_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(code_type):
            output.append((test["input"], test[code_type]))
    return output
