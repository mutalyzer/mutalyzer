from pathlib import Path

import pytest
from mutalyzer_retriever.retriever import NoReferenceError


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    try:
        with open(str(data_file), "r") as file:
            content = file.read()
    except FileNotFoundError:
        raise NoReferenceError({}, [])
    return content


def retrieve_raw(
    reference_id,
    reference_source=None,
    reference_type=None,
    size_off=True,
    configuration_path=None,
    timeout=1,
):
    if reference_type == "fasta":
        return _get_content("data/" + reference_id + ".fasta"), "fasta", "ncbi"
    elif reference_id.startswith("LRG_"):
        return _get_content("data/" + reference_id), "lrg", "lrg"
    else:
        return _get_content("data/" + reference_id + ".gff3"), "gff3", "ncbi"


def get_cds_to_mrna(cds_id, timeout=10):
    return None


@pytest.fixture(autouse=True)
def monkey_patches(monkeypatch):
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", retrieve_raw)
    monkeypatch.setattr("mutalyzer.description.get_cds_to_mrna", get_cds_to_mrna)
    monkeypatch.setattr("mutalyzer.util.configuration", lambda: None)


def code_in(code, messages):
    for meessage in messages:
        if meessage["code"] == code:
            return True
    return False


def get_codes(messages):
    return set([message["code"] for message in messages])
