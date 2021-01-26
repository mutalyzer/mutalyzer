import pytest

from normalizer.normalizer import normalize

from .commons import patch_retriever
from .test_set import TESTS_ALL


def get_tests(tests):
    output = []
    for test in tests:
        if test.get("to_test") and test["normalized"]:
            output.append((test["input"], test["normalized"]))
    return output


@pytest.mark.parametrize("input_description, normalized", get_tests(TESTS_ALL))
def test_normalizer(input_description, normalized):
    assert normalize(input_description) == normalized


def test_normalizer_other(
    i_d="LRG_303:g.[105_106del;6681G>C;6883_6884insTTTCGCCCCTTTCGCCCC]",
    n_d="LRG_303:g.[108_109del;6681G>C;6883_6884insTTTCGCCCCTTTCGCCCC]",
):
    assert normalize(i_d) == n_d
