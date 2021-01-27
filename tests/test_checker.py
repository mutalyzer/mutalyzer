import pytest

from normalizer.name_checker import name_check

from .commons import code_in

TESTS = [
    # ("NG_007485.1:g.33741del", "EOUTOFSEQBOUNDS"),
    ("NG_007485.1:g.0del", "EOUTOFBOUNDARY"),
    ("NG_007485.1:g.-1del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NM_000077.4):c.40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NM_000077.4):c.-40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NR_003529.3):n.40000del", "EOUTOFBOUNDARY"),
    ("NG_007485.1(NR_003529.3):n.-40000del", "EOUTOFBOUNDARY"),
    # ("NG_007485.1:g.4_3del", "ERANGE"),
]


@pytest.mark.parametrize("input_description, code", TESTS)
def test_messages(input_description, code):
    assert code_in(code, name_check(input_description)["errors"])
