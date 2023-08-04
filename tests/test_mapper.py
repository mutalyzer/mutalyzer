import pytest

from mutalyzer.mapper import map_description

from .commons import code_in, monkey_patches

TEST_SET = [
    (
        # Transcript to same transcript.
        ("NM_003002.2:c.274G>T", "NM_003002.2"),
        "NM_003002.2:c.274G>T",
    ),
    (
        # Transcript to same transcript.
        ("NM_003002.4:c.274G>T", "NM_003002.4"),
        "NM_003002.4:c.274G>T",
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
        (
            "NG_012337.3(NM_003002.4):c.[171del;274G>T]",
            "NG_012337.1",
            "NM_003002.2",
            "transcript",
            False,
        ),
        "NG_012337.1(NM_003002.2):c.[-60_-35del;171del;274G>T]",
    ),
    (
        (
            "NG_012337.3(NM_003002.4):c.[171del;274G>T]",
            "NG_012337.1",
            "NM_003002.2",
            "transcript",
            True,
        ),
        "NG_012337.1(NM_003002.2):c.[171del;274G>T]",
    ),
    (
        (
            "NG_012337.3(NM_003002.4):c.[171del;274G>T]",
            "NG_012337.1",
            "NM_003002.2",
            "gene",
            False,
        ),
        "NG_012337.1(NM_003002.2):c.[-60_-35del;171del;274G>T]",
    ),
    (
        (
            "NG_012337.3(NM_003002.4):c.[171del;274G>T]",
            "NG_012337.1",
            "NM_003002.2",
            "gene",
            True,
        ),
        "NG_012337.1(NM_003002.2):c.[171del;274G>T]",
    ),
    (
        # Reverse strand.
        ("NM_012459.2:c.-20del", "NG_012337.1", "NM_012459.2", "transcript", False),
        "NG_012337.1(NM_012459.2):c.-20del",
    ),
    (
        # Reverse strand.
        ("NM_012459.2:c.-20del", "NG_012337.1", "NM_012459.2", "transcript", True),
        "NG_012337.1(NM_012459.2):c.-20del",
    ),
    (
        # Reverse strand.
        ("NM_012459.2:c.-20del", "NG_012337.1", "NM_012459.2", None, False),
        "NG_012337.1(NM_012459.2):c.[-11025_-30del;-20del;129+5_133del;*496_*3448delinsA]",
    ),
    (
        # Reverse strand.
        ("NM_012459.2:c.-20del", "NG_012337.1", "NM_012459.2", "transcript", False),
        "NG_012337.1(NM_012459.2):c.-20del",
    ),
    (
        # Reverse strand.
        (
            "NG_012337.1(NM_012459.2):c.-20del",
            "NG_012337.3",
            "NM_012459.4",
            None,
            False,
        ),
        "NG_012337.3(NM_012459.4):c.[-34907_-11072del;-65del]",
    ),
    (
        # Reverse strand.
        ("NG_012337.1(NM_012459.2):c.-20del", "NG_012337.3", "NM_012459.4", None, True),
        "NG_012337.3(NM_012459.4):c.-65del",
    ),
    (
        ("NM_003002.4:c.[171del;274G>T]", "NG_012337.3", None, "gene", False),
        "NG_012337.3:g.[4_5029del;5114_6010del;6128_7022delinsC;7125G>T;7167_12959del;13949_39784del]",
    ),
    (
        ("NM_003002.4:c.[171del;274G>T]", "NG_012337.3", None, "transcript", False),
        "NG_012337.3:g.[4_5029del;5114_6010del;6128_7022delinsC;7125G>T;7167_12959del;13949_39784del]",
    ),
    (
        ("NM_003002.4:c.274G>T", "NG_012337.3", None, "transcript", True),
        "NG_012337.3:g.7125G>T",
    ),
    (
        ("NM_012459.4:c.=", "NM_012459.2", None, None, False),
        "NM_012459.2:c.[-29_13del;*496del]",
    ),
    (
        ("NM_012459.4:c.=", "NM_012459.2", None, "transcript", False),
        "NM_012459.2:c.-29_13del",
    ),
    (
        ("NM_012459.4:c.=", "NM_012459.2", "NM_012459.2", "transcript", False),
        "NM_012459.2:c.-29_13del",
    ),
    (
        ("NM_012459.2:c.=", "NM_012459.4", None, None, False),
        "NM_012459.4:c.[-33_-32insAGTCGAGAGGCGGTGCACACCCGTCGCGCATGCGCAAACACA;*495dup]",
    ),
    (
        ("NM_012459.2:c.=", "NM_012459.4", None, "transcript", False),
        "NM_012459.4:c.-33_-32insAGTCGAGAGGCGGTGCACACCCGTCGCGCATGCGCAAACACA",
    ),
    (
        ("NM_012459.2:c.=", "NM_012459.4", "NM_012459.4", "transcript", False),
        "NM_012459.4:c.-33_-32insAGTCGAGAGGCGGTGCACACCCGTCGCGCATGCGCAAACACA",
    ),
    (
        ("NG_012337.1(NM_012459.2):c.130G>A", "NM_012459.2", None, "transcript", False),
        "NM_012459.2:c.130G>A",
    ),
    (
        (
            "NG_012337.1(NM_012459.2):c.130G>A",
            "NM_012459.2",
            "NM_012459.2",
            "transcript",
            False,
        ),
        "NM_012459.2:c.130G>A",
    ),
    (
        ("NM_012459.2:c.130G>A", "NG_012337.1", "NM_012459.2", "transcript", False),
        "NG_012337.1(NM_012459.2):c.130G>A",
    ),
    (
        ("NG_012337.3(NM_003002.4):c.274G>T", "NM_003002.4", None, "transcript", False),
        "NM_003002.4:c.274G>T",
    ),
]


@pytest.mark.parametrize("input_params, correct_output", TEST_SET)
def test_mapper(input_params, correct_output):
    assert map_description(*input_params)["mapped_description"] == correct_output


TEST_ERROR = [
    (
        ("NM_003002.2:c.-31del", "NG_012337.3", "NM_003002.4", "transcript", True),
        "EMAPFILTER",
        "input",
    ),
    (
        ("NM_003002.2:c.-31del", "NG_012337.3", "NM_003002.2", "transcript", True),
        "ENOSELECTORFOUND",
        "input",
    ),
    (
        (
            "NM_003002.2:c.274G>T",
            "NG_012337.3",
            "NM_003002.4",
            "transcript",
            True,
            1000,
        ),
        "ESEQUENCELENGTH",
        "input",
    ),
    (
        ("NG_012337.3(NM_003002.4):c.274G>T", "NM_003002.4", None, None, False, 10000),
        "ESEQUENCELENGTH",
        "input",
    ),
    (
        (
            "NG_012337.3(NM_003002.4):c.[-2000del;52_53del;168_169+10del;170-2_170del;171del;274G>T]",
            "NG_012337.1",
            "NM_003002.2",
            "gene",
            False,
            10000,
        ),
        "ELOCATIONSLICE",
        "input",
    ),
    (
        (
            "NG_012337.3(NM_003002.4):c.[52_53del;168_169+10del;170-2_170del;171del;274G>T]",
            "NG_012337.1",
            "NM_003002.2",
            "transcript",
            False,
            10000,
        ),
        "ELOCATIONSLICE",
        "input",
    ),
    (
        (
            "NG_012337.3(NM_003002.4):c.[52_53del;168_169+10del;170-2_170del;171del;274G>T]",
            "NG_012337.1",
            "NM_003002.2",
            "transcr",
            False,
            10000,
        ),
        "ESLICEOPTION",
        "input",
    ),
    (
        (
            "NG_012337.3:g.7125G>T",
            "NM_003002.4",
            None,
            None,
            False,
            100000,
        ),
        "ELENGTHSDIFFERENCE",
        "input",
    ),
    (
        (
            "NG_012337.3:g.7125G>T",
            "NM_003002.4",
            None,
            "gene",
            False,
            100000,
        ),
        "ELENGTHSDIFFERENCE",
        "input",
    ),
    (
        (
            "NG_012337.3:g.7125G>T",
            "NM_003002.4",
            None,
            "transcript",
            False,
            100000,
        ),
        "ELENGTHSDIFFERENCE",
        "input",
    ),
    (
        (
            "NM_003002.4:c.169_170insA",
            "NG_012337.3",
            "NM_003002.4",
            "transcript",
            False,
            100000,
        ),
        "EINSERTIONRANGE",
        "output",
    ),
]


@pytest.mark.parametrize("input_params, error_code, source", TEST_ERROR)
def test_mapper_error(input_params, error_code, source):
    assert code_in(error_code, map_description(*input_params)["errors"])
    assert map_description(*input_params)["source"] == source
