import json
from pathlib import Path

import pytest

from normalizer.position_converter import PositionConvert

from .test_set import TESTS_ALL


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


def test_error_no_required_inputs_reference_id(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id=None,
        from_selector_id=None,
        from_coordinate_system=None,
        position=100,
        to_coordinate_system=None,
        to_selector_id=None,
        include_overlapping=False,
    )
    assert p_c.errors[0]["code"] == "ENOINPUTS"


def test_error_no_required_inputs_position(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1",
        from_selector_id=None,
        from_coordinate_system=None,
        position=None,
        to_coordinate_system=None,
        to_selector_id=None,
        include_overlapping=False,
    )
    assert p_c.errors[0]["code"] == "ENOINPUTS"


def test_error_no_required_inputs_other(monkeypatch):
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
    assert p_c.errors[0]["code"] == "ENOINPUTSOTHER"


def test_error_not_supported_reference_moltype(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NM_003002.4", from_selector_id="NM_", position=100,
    )
    assert p_c.errors[0]["code"] == "EUNSUPPORTEDREF"


def test_error_no_from_selector(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1", from_selector_id="NM_", position=100,
    )
    print(p_c.errors)
    assert p_c.errors[0]["code"] == "ENOFROMSELECTOR"


def test_info_from_selector_model_constructed_from_selector_id(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1", from_selector_id="NM_003002.2", position=100,
    )
    assert p_c.infos[0]["code"] == "IFROMSELECTOR"


def test_info_from_selector_model_constructed_from_reference(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1", position=100, to_selector_id="NM_003002.2",
    )
    assert p_c.infos[0]["code"] == "IFROMSELECTOR"


def test_info_to_selector_model_constructed_from_reference(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1",
        from_selector_id="NM_003002.2",
        from_coordinate_system="c",
        position=100,
    )
    assert p_c.infos[0]["code"] == "ITOSELECTOR"


def test_error_from_to_selectors_equal_g(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1",
        from_coordinate_system="g",
        position=100,
        to_coordinate_system="g",
    )
    assert p_c.errors[0]["code"] == "EFROMTOSELECTORSEQUAL"


def test_error_no_to_selector(monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)

    p_c = PositionConvert(
        reference_id="NG_012337.1", position=100, to_selector_id="NM_",
    )
    assert p_c.errors[0]["code"] == "ENOTOSELECTOR"
