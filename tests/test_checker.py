import json
from pathlib import Path

import pytest

from normalizer.normalizer import mutalyzer3


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), "r") as file:
        content = file.read()
    return content


def fetch_annotation(reference_id, reference_type=None):
    return _get_content("data/" + reference_id + ".gff3"), "gff3", "ncbi"


def fetch_sequence(reference_id, reference_source=None):
    return json.loads(_get_content("data/" + reference_id + ".sequence"))


TESTS = [
    ("NG_007485.1:g.33741del", "EOUTOFSEQBOUNDS"),
    ("NG_007485.1:g.0del", "EOUTOFSEQBOUNDS"),
    ("NG_007485.1:g.-1del", "EOUTOFSEQBOUNDS"),
    ("NG_007485.1(NM_000077.4):c.40000del", "EOUTOFSEQBOUNDS"),
    ("NG_007485.1(NM_000077.4):c.-40000del", "EOUTOFSEQBOUNDS"),
    ("NG_007485.1(NR_003529.3):n.40000del", "EOUTOFSEQBOUNDS"),
    ("NG_007485.1(NR_003529.3):n.-40000del", "EOUTOFSEQBOUNDS"),
    ("NG_007485.1:g.4_3del", "ERANGE"),
]


@pytest.mark.parametrize("input_description, code", TESTS)
def test_messages(input_description, code, monkeypatch):
    monkeypatch.setattr(
        "mutalyzer_retriever.retriever.fetch_annotations", fetch_annotation
    )
    monkeypatch.setattr("mutalyzer_retriever.retriever.fetch_sequence", fetch_sequence)
    messages = [list(m.keys())[0] for m in mutalyzer3(input_description)["messages"]]

    assert code in messages
