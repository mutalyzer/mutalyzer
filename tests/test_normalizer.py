from pathlib import Path

import pytest

from normalizer.normalizer import normalize

from .test_set import TESTS_ALL


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def retrieve_raw(
    reference_id,
    reference_source=None,
    reference_type=None,
    size_off=True,
    configuration_path=None,
):
    if reference_type == "fasta":
        return _get_content("data/" + reference_id + ".fasta"), "fasta", "ncbi"
    elif reference_id.startswith("LRG_"):
        return _get_content("data/" + reference_id), "lrg", "lrg"
    else:
        return _get_content("data/" + reference_id + ".gff3"), "gff3", "ncbi"


def get_tests(tests):
    output = []
    for test in tests:
        if test.get("to_test") and test["normalized"]:
            output.append((test["input"], test["normalized"]))
    return output


@pytest.mark.parametrize("input_description, normalized", get_tests(TESTS_ALL))
def test_normalizer(input_description, normalized, monkeypatch):
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", retrieve_raw)
    monkeypatch.setattr("normalizer.util.configuration", lambda: None)
    assert normalize(input_description) == normalized


def test_normalizer_other(
        monkeypatch,
        i_d='LRG_303:g.[105_106del;6681G>C;6883_6884insTTTCGCCCCTTTCGCCCC]',
        n_d='LRG_303:g.[108_109del;6681G>C;6868_6869ins[6869_6883;TTT]]',
):
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", retrieve_raw)
    monkeypatch.setattr("normalizer.util.configuration", lambda: None)
    assert normalize(i_d) == n_d
