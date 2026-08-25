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


@pytest.mark.parametrize(
    "input_description, protein_description",
    [
        ("NM_020451.3:c.100C>T", "NM_020451.3(NP_065184.2):p.?"),
        ("NM_020451.3:c.380G>A", "NM_020451.3(NP_065184.2):p.?"),
        ("NM_020451.3:c.1770del", "NM_020451.3(NP_065184.2):p.?"),
        (
            "NG_009930.1(NM_020451.3):c.100C>T",
            "NG_009930.1(NP_065184.2):p.?",
        ),
    ],
)
def test_normalize_alt_in_frame_stop_codon(input_description, protein_description):
    d = normalize_alt(input_description)
    assert [info["code"] for info in d["infos"]] == ["IINFRAMESTOPCODON"]
    assert d["protein"]["description"] == protein_description
