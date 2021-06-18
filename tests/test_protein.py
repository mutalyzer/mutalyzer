import pytest

from normalizer.name_checker import name_check

from .commons import code_in, patch_retriever

TESTS = [
    {
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.(=)",
        "normalized": "NG_012337.1(NP_002993.1):p.(=)",
        "to_test": True,
    },
    {
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.=",
        "normalized": "NG_012337.1(NP_002993.1):p.=",
        "to_test": True,
    },
    {
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92=",
        "normalized": "NG_012337.1(NP_002993.1):p.Asp92=",
        "to_test": True,
    },
    {
        "keywords": ["protein", "substitution", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92Tyr",
        "normalized": "NG_012337.1(NP_002993.1):p.(Asp92Tyr)",
        "to_test": True,
    },
    {
        "keywords": ["protein", "delins", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92delinsTyr",
        "normalized": "NG_012337.1(NP_002993.1):p.(Asp92Tyr)",
        "to_test": True,
    },
]


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test"):
            if test.get(t_type):
                output.append((test["input"], test[t_type]))
            else:
                output.append((test["input"], None))
    return output


@pytest.mark.parametrize(
    "input_description, normalized",
    get_tests(TESTS, "normalized"),
)
def test_protein(input_description, normalized):
    d = name_check(input_description)
    assert d["normalized_description"] == normalized
