import pytest
from mutalyzer_hgvs_parser import to_model

from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer.converter.to_internal_coordinates import to_internal_coordinates
from mutalyzer.converter.to_internal_indexing import to_internal_indexing
from mutalyzer.description_model import model_to_string

from .generator import append_transcript, generate_references

TESTS_SET = [
    {
        "c": {
            "id": "t1",
            "type": "mRNA",
            "inverted": False,
            "exon": [(3, 6), (8, 13), (16, 21), (24, 26)],
            "cds": (10, 19),
        },
        "n": {
            "id": "t2",
            "type": "ncRNA",
            "inverted": False,
            "exon": [(3, 6), (8, 13), (16, 21), (24, 26)],
        },
        "deleted": [
            {
                "i": "0_1",
                "g": {"hgvs": "1"},
                "c": {"hgvs": "-8", "other": ["-5-3"]},
                "n": {"other": ["1-3", "-3"]},
            },
            {
                "i": "1_2",
                "g": {"hgvs": "2"},
                "c": {"hgvs": "-7", "other": ["-5-2"]},
                "n": {"hgvs": "1-2", "other": ["-2"]},
            },
            {
                "i": "2_3",
                "g": {"hgvs": "3"},
                "c": {"hgvs": "-6", "other": ["-5-1"]},
                "n": {"hgvs": "1-1", "other": ["-1"]},
            },
            {
                "i": "3_4",
                "g": {"hgvs": "4"},
                "c": {"hgvs": "-5", "other": []},
                "n": {"hgvs": "1", "other": []},
            },
            {
                "i": "4_5",
                "g": {"hgvs": "5"},
                "c": {"hgvs": "-4", "other": []},
                "n": {"hgvs": "2", "other": []},
            },
            {
                "i": "5_6",
                "g": {"hgvs": "6"},
                "c": {"hgvs": "-3", "other": []},
                "n": {"hgvs": "3", "other": []},
            },
            {
                "i": "6_7",
                "g": {"hgvs": "7"},
                "c": {"hgvs": "-3+1", "other": ["-2-2"]},
                "n": {"hgvs": "3+1", "other": ["4-2"]},
            },
            {
                "i": "7_8",
                "g": {"hgvs": "8"},
                "c": {"hgvs": "-2-1", "other": []},
                "n": {"hgvs": "4-1", "other": []},
            },
            {
                "i": "8_9",
                "g": {"hgvs": "9"},
                "c": {"hgvs": "-2", "other": []},
                "n": {"hgvs": "4", "other": []},
            },
            {
                "i": "9_10",
                "g": {"hgvs": "10"},
                "c": {"hgvs": "-1", "other": []},
                "n": {"hgvs": "5", "other": []},
            },
            {
                "i": "10_11",
                "g": {"hgvs": "11"},
                "c": {"hgvs": "1", "other": ["1"]},
                "n": {"hgvs": "6", "other": ["6"]},
            },
            {
                "i": "11_12",
                "g": {"hgvs": "12"},
                "c": {"hgvs": "2", "other": ["2"]},
                "n": {"hgvs": "7", "other": ["7"]},
            },
            {
                "i": "12_13",
                "g": {"hgvs": "13"},
                "c": {"hgvs": "3", "other": ["3"]},
                "n": {"hgvs": "8", "other": ["8"]},
            },
            {
                "i": "13_14",
                "g": {"hgvs": "14"},
                "c": {"hgvs": "3+1", "other": ["4-3"]},
                "n": {"hgvs": "8+1", "other": ["9-3"]},
            },
            {
                "i": "14_15",
                "g": {"hgvs": "15"},
                "c": {"hgvs": "3+2", "other": ["4-2"]},
                "n": {"hgvs": "8+2", "other": ["9-2"]},
            },
            {
                "i": "15_16",
                "g": {"hgvs": "16"},
                "c": {"hgvs": "4-1", "other": ["3+3"]},
                "n": {"hgvs": "9-1", "other": ["8+3"]},
            },
            {
                "i": "16_17",
                "g": {"hgvs": "17"},
                "c": {"hgvs": "4", "other": ["4"]},
                "n": {"hgvs": "9", "other": ["9"]},
            },
            {
                "i": "17_18",
                "g": {"hgvs": "18"},
                "c": {"hgvs": "5", "other": ["5"]},
                "n": {"hgvs": "10", "other": ["10"]},
            },
            {
                "i": "18_19",
                "g": {"hgvs": "19"},
                "c": {"hgvs": "6", "other": ["6"]},
                "n": {"hgvs": "11", "other": ["11"]},
            },
            {
                "i": "19_20",
                "g": {"hgvs": "20"},
                "c": {"hgvs": "*1", "other": ["*1"]},
                "n": {"hgvs": "12", "other": ["12"]},
            },
            {
                "i": "20_21",
                "g": {"hgvs": "21"},
                "c": {"hgvs": "*2", "other": ["*2"]},
                "n": {"hgvs": "13", "other": ["13"]},
            },
            {
                "i": "21_22",
                "g": {"hgvs": "22"},
                "c": {"hgvs": "*2+1", "other": ["*2+1"]},
                "n": {"hgvs": "13+1", "other": ["13+1"]},
            },
            {
                "i": "22_23",
                "g": {"hgvs": "23"},
                "c": {"hgvs": "*2+2", "other": ["*2+2"]},
                "n": {"hgvs": "13+2", "other": ["13+2"]},
            },
            {
                "i": "23_24",
                "g": {"hgvs": "24"},
                "c": {"hgvs": "*3-1", "other": ["*2+3"]},
                "n": {"hgvs": "14-1", "other": ["13+3"]},
            },
            {
                "i": "24_25",
                "g": {"hgvs": "25"},
                "c": {"hgvs": "*3", "other": ["*3"]},
                "n": {"hgvs": "14", "other": ["14"]},
            },
            {
                "i": "25_26",
                "g": {"hgvs": "26"},
                "c": {"hgvs": "*4", "other": ["*4"]},
                "n": {"hgvs": "15", "other": ["15"]},
            },
            {
                "i": "26_27",
                "g": {"hgvs": "27"},
                "c": {"hgvs": "*5", "other": ["*4+1"]},
                "n": {"other": ["15+1", "*1"]},
            },
            {
                "i": "27_28",
                "g": {"hgvs": "28"},
                "c": {"hgvs": "*6", "other": ["*4+2"]},
                "n": {"other": ["15+2", "*2"]},
            },
            {
                "i": "28_29",
                "g": {"hgvs": "29"},
                "c": {"hgvs": "*7", "other": ["*4+3"]},
                "n": {"other": ["15+3", "*3"]},
            },
            {
                "i": "3_6",
                "g": {"hgvs": "4_6"},
                "c": {"hgvs": "-5_-3", "other": []},
                "n": {"hgvs": "1_3", "other": []},
            },
            {
                "i": "3_7",
                "g": {"hgvs": "4_7"},
                "c": {"hgvs": "-5_-3+1", "other": []},
                "n": {"hgvs": "1_3+1", "other": []},
            },
            {
                "i": "2_7",
                "g": {"hgvs": "3_7"},
                "c": {"hgvs": "-6_-3+1", "other": []},
                "n": {"hgvs": "1-1_3+1", "other": []},
            },
            {
                "i": "7_14",
                "g": {"hgvs": "8_14"},
                "c": {"hgvs": "-2-1_3+1", "other": []},
                "n": {"hgvs": "4-1_8+1", "other": []},
            },
            {
                "i": "7_22",
                "g": {"hgvs": "8_22"},
                "c": {"hgvs": "-2-1_*2+1", "other": []},
                "n": {"hgvs": "4-1_13+1", "other": []},
            },
            {
                "i": "0_29",
                "g": {"hgvs": "1_29"},
                "c": {"hgvs": "-8_*7", "other": []},
                "n": {"hgvs": "1-3_15+3", "other": []},
            },
            {
                "i": "0_?",
                "g": {"hgvs": "1_?"},
                "c": {"hgvs": "-8_?", "other": []},
                "n": {"hgvs": "1-3_?", "other": []},
            },
            {
                "i": "(0_?)",
                "g": {"hgvs": "(1_?)"},
                "c": {"hgvs": "(-8_?)", "other": []},
                "n": {"hgvs": "(1-3_?)", "other": []},
            },
            {
                "i": "(0_2)",
                "g": {"hgvs": "(1_2)"},
                "c": {"hgvs": "(-8_-7)", "other": []},
                "n": {"hgvs": "(1-3_1-2)", "other": []},
            },
        ],
        "deleted_insertion": [
            {
                "i": "1_1",
                "g": {"hgvs": "1_2"},
                "c": {"hgvs": "-8_-7", "other": ["-7-1_-5-2"]},
                "n": {"hgvs": "1-3_1-2", "other": ["*3-28_13-19"]},
            },
        ],
    },
    {
        "c": {
            "id": "t1",
            "type": "mRNA",
            "inverted": False,
            "exon": [(3, 4)],
            "cds": (3, 4),
        },
        "n": {"id": "t2", "type": "ncRNA", "inverted": False, "exon": [(3, 4)]},
        "deleted": [
            {
                "i": "0_1",
                "g": {"hgvs": "1"},
                "c": {"hgvs": "-3", "other": ["-2-1", "*1-4"]},
                "n": {"hgvs": "-3", "other": ["-3"]},
            },
            {
                "i": "2_3",
                "g": {"hgvs": "3"},
                "c": {"hgvs": "-1", "other": ["-1", "*1-2"]},
                "n": {"hgvs": "-1", "other": ["-1"]},
            },
            {
                "i": "3_4",
                "g": {"hgvs": "4"},
                "c": {"hgvs": "1", "other": ["2-1", "-1+1"]},
                "n": {"hgvs": "1", "other": ["1"]},
            },
        ],
    },
]


def generate_to_internal_location_test(t, loc, d, refs):
    tests = [
        (
            d.format("", "g", loc["g"]["hgvs"]),
            d.format("", "i", loc["i"]),
            refs,
        ),
        (
            d.format("({})".format(t["c"]["id"]), "c", loc["c"]["hgvs"]),
            d.format("", "i", loc["i"]),
            refs,
        ),
    ]
    if loc["c"].get("other"):
        for other in loc["c"]["other"]:
            tests.append(
                (
                    d.format("({})".format(t["c"]["id"]), "c", other),
                    d.format("", "i", loc["i"]),
                    refs,
                )
            )
    if loc["n"].get("hgvs"):
        tests.append(
            (
                d.format("({})".format(t["n"]["id"]), "n", loc["n"]["hgvs"]),
                d.format("", "i", loc["i"]),
                refs,
            ),
        )
    if loc["n"].get("other"):
        for other in loc["n"]["other"]:
            tests.append(
                (
                    d.format("({})".format(t["n"]["id"]), "n", other),
                    d.format("", "i", loc["i"]),
                    refs,
                )
            )
    return tests


def generate_to_internal_locations_tests(tests_set):
    deleted = [
        "R1{}:{}.{}A>T",
        "R1{}:{}.{}del",
        "R1{}:{}.{}delinsA",
        "R1{}:{}.{}dup",
        "R1{}:{}.{}inv",
        "R1{}:{}.{}=",
    ]
    deleted_insertion = [
        "R1{}:{}.{}insA",
    ]

    tests = []
    for test in tests_set:
        references = generate_references(test["c"])
        append_transcript(references, test["n"])
        for description in deleted:
            for location in test["deleted"]:
                tests += generate_to_internal_location_test(
                    test, location, description, references
                )
        for description in deleted_insertion:
            for location in test["deleted_insertion"]:
                tests += generate_to_internal_location_test(
                    test, location, description, references
                )
        return tests


@pytest.mark.parametrize(
    "description_in, description_expected, references",
    generate_to_internal_locations_tests(TESTS_SET),
)
def test_to_internal_locations(description_in, description_expected, references):
    description_model = to_model(description_in)
    description_out = model_to_string(
        to_internal_indexing(to_internal_coordinates(description_model, references))
    )
    assert description_out == description_expected


def generate_to_hgvs_location_test(t, loc, d, refs):
    tests = [
        (
            d.format("", "i", loc["i"]),
            d.format("", "g", loc["g"]["hgvs"]),
            refs,
        ),
        (
            d.format("", "i", loc["i"]),
            d.format("({})".format(t["c"]["id"]), "c", loc["c"]["hgvs"]),
            refs,
        ),
    ]
    if loc["n"].get("hgvs"):
        tests.append(
            (
                d.format("", "i", loc["i"]),
                d.format("({})".format(t["n"]["id"]), "n", loc["n"]["hgvs"]),
                refs,
            )
        )

    return tests


def generate_to_hgvs_locations_tests(tests_set):
    deleted = [
        "R1{}:{}.{}A>T",
        "R1{}:{}.{}del",
        "R1{}:{}.{}delinsA",
        "R1{}:{}.{}dup",
        "R1{}:{}.{}inv",
        "R1{}:{}.{}=",
    ]
    deleted_insertion = [
        "R1{}:{}.{}insA",
    ]

    tests = []
    for test in tests_set:
        references = generate_references(test["c"])
        append_transcript(references, test["n"])
        for description in deleted:
            for location in test["deleted"]:
                tests += generate_to_hgvs_location_test(
                    test, location, description, references
                )
        for description in deleted_insertion:
            for location in test["deleted_insertion"]:
                tests += generate_to_hgvs_location_test(
                    test, location, description, references
                )
        return tests


@pytest.mark.parametrize(
    "description_in, description_expected, references",
    generate_to_hgvs_locations_tests(TESTS_SET),
)
def test_to_hgvs_locations(description_in, description_expected, references):
    model_internal_indexing = to_model(description_in)
    model_hgvs = to_model(description_expected)
    if model_hgvs["coordinate_system"] in ["c", "n"]:
        to_selector_id = model_hgvs["reference"]["selector"]["id"]
    else:
        to_selector_id = None

    hgvs_conversion = model_to_string(
        to_hgvs_locations(
            model=model_internal_indexing,
            references=references,
            to_coordinate_system=model_hgvs["coordinate_system"],
            to_selector_id=to_selector_id,
            degenerate=True,
        )
    )

    assert hgvs_conversion == description_expected


EQUIVALENT_DESCRIPTIONS = [
    # Deletion
    ("R1:g.4del", "R1:i.3_4del"),
    ("R1:g.3_4del", "R1:i.2_4del"),
    # Substitution
    ("R1:g.4A>T", "R1:i.3_4A>T"),
    ("R1:g.4_6AAT>G", "R1:i.3_6AAT>G"),
    ("R1:g.4_6AAT>6_9", "R1:i.3_6AAT>5_9"),
    ("R1(t1):c.1A>T", "R1:i.10_11A>T"),
    ("R1(t1):c.4A>T", "R1:i.16_17A>T"),
    ("R1(t1):c.2_6AAT>G", "R1:i.11_19AAT>G"),
    ("R1(t1):c.4_6AAT>2_4", "R1:i.16_19AAT>11_17"),
    # Insertion
    ("R1:g.4_5insT", "R1:i.4_4insT"),
    ("R1:g.4_5ins7_8", "R1:i.4_4ins6_8"),
    ("R1:g.4_5ins[7_8;10_20]", "R1:i.4_4ins[6_8;9_20]"),
    #
    # Duplication
    ("R1:g.4dup", "R1:i.3_4dup"),
    ("R1:g.3_4dup", "R1:i.2_4dup"),
    # # Inversion
    ("R1:g.11_13inv", "R1:i.10_13inv"),
    # Conversion
    ("R1:g.4_5con7_8", "R1:i.3_5con6_8"),
    #
    # Deletion-insertion
    ("R1:g.4delins7_8", "R1:i.3_4delins6_8"),
    ("R1:g.4delins6_31", "R1:i.3_4delins5_31"),
    ("R1(t1):c.-8_*7A>T", "R1:i.0_29A>T"),
    ("R1:g.(?_?)del", "R1:i.(?_?)del"),
    ("R1:g.?_?del", "R1:i.?_?del"),
    ("R1:g.(?_?)_(?_?)del", "R1:i.(?_?)_(?_?)del"),
    ("R1:g.3_?del", "R1:i.2_?del"),
    ("R1:g.(4_7)del", "R1:i.(3_7)del"),
    ("R1:g.(4_?)del", "R1:i.(3_?)del"),
    ("R1:g.(4_?)_(5_7)del", "R1:i.(3_?)_(4_7)del"),
    ("R1:g.(4_?)_7del", "R1:i.(3_?)_7del"),
    ("R1:g.(?_4)_7del", "R1:i.(?_4)_7del"),
    ("R1:g.4_(5_?)del", "R1:i.3_(4_?)del"),
    ("R1(t1):c.-4A>T", "R1:i.4_5A>T"),
    ("R1(t1):c.-2-1A>T", "R1:i.7_8A>T"),
    ("R1(t1):c.-2-1_*2+1del", "R1:i.7_22del"),
    ("R1(t2):n.1-3_15+3A>T", "R1:i.0_29A>T"),
]

TO_INTERNAL_ONLY = [
    ("R1:g.?del", "R1:i.?_?del"),
    ("R1(t2):n.*1A>T", "R1:i.26_27A>T"),
    ("R1(t2):n.-1A>T", "R1:i.2_3A>T"),
    ("R1(t2):n.-2A>T", "R1:i.1_2A>T"),
    ("R1(t2):n.-3A>T", "R1:i.0_1A>T"),
    ("R1:g.4_5ins[R1:g.7_8;10_20]", "R1:i.4_4ins[R1:i.6_8;9_20]"),
    ("R1:g.4_5ins[R1(t1):c.7_8;10_20]", "R1:i.4_4ins[R1:i.19_21;9_20]"),
]


@pytest.mark.parametrize(
    "hgvs, hgvs_internal_indexing",
    EQUIVALENT_DESCRIPTIONS + TO_INTERNAL_ONLY,
)
def test_to_internal_coordinates_simple(hgvs, hgvs_internal_indexing):
    d_m = to_model(hgvs)
    r_model = generate_references(
        {
            "id": "t1",
            "type": "mRNA",
            "inverted": False,
            "exon": [(3, 6), (8, 13), (16, 21), (24, 26)],
            "cds": (10, 19),
        }
    )
    append_transcript(
        r_model,
        {
            "id": "t2",
            "type": "ncRNA",
            "inverted": False,
            "exon": [(3, 6), (8, 13), (16, 21), (24, 26)],
        },
    )

    hgvs_internal_indexing_conversion = model_to_string(
        to_internal_indexing(to_internal_coordinates(d_m, r_model))
    )

    assert hgvs_internal_indexing_conversion == hgvs_internal_indexing


@pytest.mark.parametrize(
    "hgvs, hgvs_internal_indexing",
    EQUIVALENT_DESCRIPTIONS,
)
def test_to_hgvs_locations_simple(hgvs, hgvs_internal_indexing):
    r_model = generate_references(
        {
            "id": "t1",
            "type": "mRNA",
            "inverted": False,
            "exon": [(3, 6), (8, 13), (16, 21), (24, 26)],
            "cds": (10, 19),
        }
    )
    append_transcript(
        r_model,
        {
            "id": "t2",
            "type": "ncRNA",
            "inverted": False,
            "exon": [(3, 6), (8, 13), (16, 21), (24, 26)],
        },
    )
    model_internal_indexing = to_model(hgvs_internal_indexing)
    model_hgvs = to_model(hgvs)
    if model_hgvs["coordinate_system"] in ["c", "n"]:
        to_selector_id = model_hgvs["reference"]["selector"]["id"]
    else:
        to_selector_id = None

    hgvs_conversion = model_to_string(
        to_hgvs_locations(
            model=model_internal_indexing,
            references=r_model,
            to_coordinate_system=model_hgvs["coordinate_system"],
            to_selector_id=to_selector_id,
            degenerate=True,
        )
    )

    assert hgvs_conversion == hgvs
