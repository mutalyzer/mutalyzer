from pathlib import Path

import pytest


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
    timeout=1,
):
    if reference_type == "fasta":
        return _get_content("data/" + reference_id + ".fasta"), "fasta", "ncbi"
    elif reference_id.startswith("LRG_"):
        return _get_content("data/" + reference_id), "lrg", "lrg"
    else:
        return _get_content("data/" + reference_id + ".gff3"), "gff3", "ncbi"


@pytest.fixture(autouse=True)
def patch_retriever(monkeypatch):
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", retrieve_raw)
    monkeypatch.setattr("normalizer.util.configuration", lambda: None)


def code_in_errors(code, errors):
    for error in errors:
        if error["code"] == code:
            return True
    return False


def get_errors_codes(errors):
    return set([error["code"] for error in errors])
