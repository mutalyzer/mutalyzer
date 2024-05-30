import pytest

from mutalyzer.spdi_converter import spdi_converter

from .commons import code_in, monkey_patches
from .variants_set import TESTS_ALL


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test") and all([test.get(k, False) for k in ["input_spdi", t_type]]):
            for spdi in test["input_spdi"]:
                output.append((spdi, test[t_type]))
    return output


@pytest.mark.parametrize(
    "input_description, normalized",
    get_tests(TESTS_ALL, "normalized")
    + [
        ("NG_012337.3(NM_003002.4):273:T:A", "NG_012337.3:g.274T>A"),
        ("NG_012337.3(NM_003002.4):273:1:A", "NG_012337.3:g.274T>A"),
        ("NM_003002.4:308:G:T", "NM_003002.4:c.274G>T"),
        ("NM_003002.4:308:1:T", "NM_003002.4:c.274G>T"),
        ("NR_038420.1:10::T", "NR_038420.1:n.10_11insT"),
    ]
)
def test_spdi_converter(input_description, normalized):
    d = spdi_converter(input_description)
    assert d["normalized_description"] == normalized


@pytest.mark.parametrize(
    "input_description, codes",
    get_tests(TESTS_ALL, "errors")
    + [
        ("NG_012337.3:71250:1:", ["EOUTOFBOUNDARY", "EOUTOFBOUNDARY"]),
        ("NG_012337.3:1:71250:", ["EOUTOFBOUNDARY"]),
        ("NG_012337.3:7125:100000:", ["EOUTOFBOUNDARY"]),
        ("NG_012337.3:7125:", ["ESYNTAXUEOF"]),
        ("NG_012337.3:7125:$", ["ESYNTAXUC"]),
    ],
)
def test_errors(input_description, codes):
    assert codes == [error["code"] for error in spdi_converter(input_description)["errors"]]
