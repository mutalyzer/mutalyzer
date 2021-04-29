import pytest

from normalizer.name_checker import name_check
from normalizer.reference import retrieve_reference, slice_to_selector

from .commons import code_in, patch_retriever

TESTS = [
    {
        "keywords": ["rna", "equal", "genomic"],
        "input": "NG_012337.1(NM_003002.2):r.275=",
        "normalized": "NG_012337.1(NM_003002.2):r.=",
        "to_test": True,
    },
    {
        "keywords": ["rna", "uncertain", "genomic"],
        "input": "NG_012337.1(NM_003002.2):r.?",
        "normalized": "NG_012337.1(NM_003002.2):r.?",
        "to_test": False,
    },
    {
        "keywords": ["rna", "substitution", "genomic"],
        "input": "NG_012337.1(NM_003002.2):r.275a>c",
        "normalized": "NG_012337.1(NM_003002.2):r.275a>c",
        "to_test": True,
    },
    {
        "keywords": ["rna", "substitution", "mRNA"],
        "input": "NM_003002.2:r.275a>c",
        "normalized": "NM_003002.2:r.275a>c",
        "to_test": True,
    },
    {
        "keywords": ["rna", "substitution", "ncRNA"],
        "input": "NR_038420.1:r.75g>a",
        "normalized": "NR_038420.1:r.75g>a",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "genomic"],
        "input": "NG_012337.1(NM_003002.2):r.273del",
        "normalized": "NG_012337.1(NM_003002.2):r.274del",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "mRNA"],
        "input": "NM_003002.2:r.273del",
        "normalized": "NM_003002.2:r.274del",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "ncRNA"],
        "input": "NR_038420.1:r.74del",
        "normalized": "NR_038420.1:r.75del",
        "to_test": True,
    },
]


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(t_type):
            output.append((test["input"], test[t_type]))
    return output


@pytest.mark.parametrize(
    "input_description, normalized",
    get_tests(TESTS, "normalized"),
)
def test_rna(input_description, normalized):
    d = name_check(input_description)
    assert d["normalized_description"] == normalized
