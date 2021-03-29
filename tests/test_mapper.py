import pytest

from normalizer.mapper import map_description

from .commons import code_in, patch_retriever


TEST_SET = [
    (
        # Transcript to same transcript.
        ("NM_003002.2:c.274G>T", "NM_003002.2"),
        "NM_003002.2:c.274G>T",
    ),
    (
        # Transcript to NG with transcript slice.
        ("NM_003002.2:c.274G>T", "NG_012337.1", "NM_003002.2", "transcript", False),
        "NG_012337.1(NM_003002.2):c.[274G>T;*824A[18]]",
    ),
    (
        # Transcript to NG with transcript slice with filtering.
        ("NM_003002.2:c.274G>T", "NG_012337.1", "NM_003002.2", "transcript", True),
        "NG_012337.1(NM_003002.2):c.274G>T",
    ),
    (
        # Transcript to NG with transcript slice updated version.
        ("NM_003002.2:c.-31del", "NG_012337.3", "NM_003002.4", "transcript", False),
        "NG_012337.3(NM_003002.4):c.[-34_-32delinsTGGGAATTGTCGCCTAAGTGGTTCCGGG;*824A[18]]",
    ),
    (
        # Reverse strand.
        ("NM_012459.2:c.-20del", "NG_012337.1", "NM_012459.2", "transcript", False),
        "NG_012337.1(NM_012459.2):c.[-20del;*496dup]",
    ),
    (
        # Reverse strand.
        ("NM_012459.2:c.-20del", "NG_012337.1", "NM_012459.2", "transcript", True),
        "NG_012337.1(NM_012459.2):c.-20del",
    ),
]


@pytest.mark.parametrize("input_params, correct_output", TEST_SET)
def test_mapper(input_params, correct_output):
    print(input_params)
    print(correct_output)
    assert map_description(*input_params) == correct_output


TEST_ERROR = [
    (
        ("NM_003002.2:c.-31del", "NG_012337.3", "NM_003002.4", "transcript", True),
        "EMAPFILTER",
    ),
    (
        ("NM_003002.2:c.-31del", "NG_012337.3", "NM_003002.2", "transcript",
         True),
        "ENOSELECTORFOUND",
    ),
]


@pytest.mark.parametrize("input_params, error_code", TEST_ERROR)
def test_mapper_error(input_params,error_code):
    print(input_params)
    print(error_code)
    assert code_in(error_code, map_description(*input_params)["errors"])
