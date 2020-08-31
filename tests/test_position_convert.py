import pytest
from .test_set import TESTS_ALL
from normalizer.position_convert import PositionConvert
from pathlib import Path
import json


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def fetch_annotation(reference_id, reference_type=None):
    return _get_content("data/" + reference_id + ".gff3"), "gff3", "ncbi"


def fetch_sequence(reference_id, reference_source=None):
    return json.loads(_get_content("data/" + reference_id + ".sequence"))


def get_tests(tests):
    output = []
    for test in tests:
        if test.get("to_test") and test["normalized"]:
            output.append((test["input"], test["normalized"]))
    return output


def test_error_no_inputs(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1",
        from_selector_id=None,
        from_coordinate_system=None,
        position=100,
        to_coordinate_system=None,
        to_selector_id=None,
        include_overlapping=False,
    )
    assert p_c.errors == [{"code": "ENOFROMSELECTOR"}]


def test_error_no_from_selector(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1",
        from_selector_id="NM_",
        from_coordinate_system=None,
        position=100,
        to_coordinate_system=None,
        to_selector_id=None,
        include_overlapping=False,
    )
    assert p_c.errors == [{"code": "ENOFROMSELECTOR"}]


def test_error_no_to_selector(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1",
        from_selector_id=None,
        from_coordinate_system=None,
        position=100,
        to_coordinate_system=None,
        to_selector_id="NM_",
        include_overlapping=False,
    )
    assert p_c.errors == [{"code": "ENOTOSELECTOR"}]
