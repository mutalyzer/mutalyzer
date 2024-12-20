import pytest

from mutalyzer.normalizer import normalize

from .commons import code_in, monkey_patches

TESTS = [
    {
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.=",
        "normalized": "NG_012337.1(NP_002993.1):p.=",
        "to_test": True,
    },
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
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92Asp",
        "normalized": "NG_012337.1(NP_002993.1):p.=",
        "to_test": True,
    },
    {
        "keywords": ["protein", "substitution", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92Tyr",
        "normalized": "NG_012337.1(NP_002993.1):p.Asp92Tyr",
        "to_test": True,
    },
    {
        "keywords": ["protein", "delins", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92delinsTyr",
        "normalized": "NG_012337.1(NP_002993.1):p.Asp92Tyr",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
            "delins",
            "substitution",
            "mRNA",
            "translation exception",
        ],
        "input": "NM_005410.4:p.U59delinsTyr",
        "normalized": "NM_005410.4(NP_005401.3):p.Sec59Tyr",
        "to_test": True,
    },
    {
        "keywords": ["protein", "delins", "mRNA", "translation exception", "error"],
        "input": "NM_005410.4:p.D59delinsA",
        "errors": [
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "D not found in the reference sequence, found U instead.",
                "paths": [["variants", 0, "location"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein", "delins", "protein", "translation exception", "error"],
        "input": "NP_005401.3:p.D59delinsA",
        "errors": [
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "D not found in the reference sequence, found U instead.",
                "paths": [["variants", 0, "location"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
            "delins",
            "mRNA",
            "reference insertion",
            "translation exception",
        ],
        "input": "NM_005410.4:p.U59delinsNP_002993.1:H50_53",
        "normalized": "NM_005410.4(NP_005401.3):p.Sec59delinsHisLeuSerPro",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
            "delins",
            "mRNA",
            "reference insertion",
            "translation exception",
            "error",
        ],
        "input": "NM_005410.4:p.D59delinsNP_002993.1:Asp50_53",
        "errors": [
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "D not found in the reference sequence, found U instead.",
                "paths": [["variants", 0, "location"]],
            },
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "D not found in the reference sequence, found H instead.",
                "paths": [["variants", 0, "inserted", 0, "location", "start"]],
            },
        ],
        "to_test": True,
    },
    {
        "keywords": [
            "protein missing mrna for backtranslation",
        ],
        "input": "YP_009725300.1:p.(Leu360Ter)",
        "normalized": "YP_009725300.1:p.(Leu360Ter)",
        "to_test": True,
    },
    {
        "keywords": [
            "protein reverse strand",
        ],
        "input": "NG_008835.1(NP_071436.1):p.(Lys34Val)",
        "normalized": "NG_008835.1(NP_071436.1):p.(Lys34Val)",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
        ],
        "input": "NP_071436.1:p.(Lys34Val)",
        "normalized": "NP_071436.1:p.(Lys34Val)",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
        ],
        "input": "NP_071436.1:p.(Lys34Xaa)",
        "normalized": "NP_071436.1:p.(Lys34Xaa)",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
        ],
        "input": "NP_071436.1:p.(K34X)",
        "normalized": "NP_071436.1:p.(Lys34Xaa)",
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
    d = normalize(input_description)
    assert d.get("normalized_description") == normalized


@pytest.mark.parametrize(
    "input_description, errors",
    get_tests(TESTS, "errors"),
)
def test_protein_errors(input_description, errors):
    d = normalize(input_description)
    assert d.get("errors") == errors
