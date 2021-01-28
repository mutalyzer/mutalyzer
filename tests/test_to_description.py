import pytest
from mutalyzer_hgvs_parser import parse_description_to_model

from normalizer.description_model import model_to_string

from .variants_set import TESTS_ALL

TESTS = [
    (
        {
            "reference": {"id": "R1"},
            "coordinate_system": "g",
            "variants": [
                {
                    "location": {
                        "type": "range",
                        "start": {
                            "type": "range",
                            "start": {"type": "point", "position": 3},
                            "end": {"type": "point", "uncertain": True},
                            "uncertain": True,
                        },
                        "end": {"type": "point", "position": 7},
                    },
                    "type": "deletion",
                    "source": "reference",
                }
            ],
        },
        "R1:g.(3_?)_7del",
    ),
    (
        {
            "reference": {"id": "R1", "selector": {"id": "S1"}},
            "coordinate_system": "c",
            "variants": [
                {
                    "location": {
                        "type": "range",
                        "start": {
                            "type": "range",
                            "start": {"type": "point", "position": 3},
                            "end": {"type": "point", "uncertain": True},
                            "uncertain": True,
                        },
                        "end": {"type": "point", "position": 7},
                    },
                    "type": "deletion",
                    "source": "reference",
                }
            ],
        },
        "R1(S1):c.(3_?)_7del",
    ),
    (
        {
            "reference": {"id": "R1", "selector": {"id": "S1"}},
            "coordinate_system": "c",
            "variants": [
                {
                    "location": {
                        "type": "range",
                        "start": {
                            "type": "range",
                            "start": {"type": "point", "position": 3},
                            "end": {"type": "point", "uncertain": True},
                            "uncertain": True,
                        },
                        "end": {
                            "type": "range",
                            "start": {"type": "point", "uncertain": True},
                            "end": {
                                "type": "point",
                                "outside_cds": "upstream",
                                "position": 10,
                            },
                            "uncertain": True,
                        },
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                    "deleted": [{"sequence": "TTAA", "source": "description"}],
                    "inserted": [
                        {
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": 59},
                                "end": {"type": "point", "position": 79},
                            },
                            "source": "reference",
                        },
                        {
                            "location": {
                                "type": "range",
                                "start": {
                                    "type": "point",
                                    "outside_cds": "downstream",
                                    "position": 3,
                                },
                                "end": {
                                    "type": "point",
                                    "position": 100,
                                    "offset": {"value": 40},
                                },
                            },
                            "source": "reference",
                        },
                    ],
                }
            ],
        },
        "R1(S1):c.(3_?)_(?_-10)delTTAAins[59_79;*3_100+40]",
    ),
]


@pytest.mark.parametrize("model, description", TESTS)
def test_model_to_string_simple(model, description):
    assert model_to_string(model) == description


def extract_descriptions():
    """
    Extract descriptions from variants_set.TEST_ALL.
    """
    variants = []
    for test in TESTS_ALL:
        for k in ["input", "normalized", "genomic"]:
            if test.get(k):
                variants.append(test[k])
    return variants


@pytest.mark.parametrize(
    "description",
    extract_descriptions(),
)
def test__model_to_string_descriptions_test_all(description):
    """
    Check if the parsed descriptions from variants_set.TEST_ALL are converted back.
    """
    assert model_to_string(parse_description_to_model(description)) == description
