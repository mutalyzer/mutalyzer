import pytest
from .test_set import TESTS_ALL
from normalizer.normalizer import mutalyzer3
from pathlib import Path
import json


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), 'r') as file:
        content = file.read()
    return content


def fetch_annotation(reference_id, reference_type=None):
    return _get_content('data/' + reference_id + '.gff3'), 'gff3', 'ncbi'


def fetch_sequence(reference_id, reference_source=None):
    return json.loads(_get_content('data/' + reference_id + '.sequence'))


def get_tests(tests):
    output = []
    for test in tests:
        if test.get('to_test') and test['normalized']:
            output.append((test['input'], test['normalized']))
    return output


@pytest.mark.parametrize('input_description, normalized',
                         get_tests(TESTS_ALL))
def test_mutalyzer3(input_description, normalized, monkeypatch):
    monkeypatch.setattr('retriever.retriever.fetch_annotations',
                        fetch_annotation)
    monkeypatch.setattr('retriever.retriever.fetch_sequence',
                        fetch_sequence)
    assert mutalyzer3(input_description)['normalized description'] == normalized
