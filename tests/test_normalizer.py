import pytest
from normalizer.normalizer import mutalyzer3
from pathlib import Path


def _get_content(relative_location):
    data_file = Path(__file__).parent.joinpath(relative_location)
    with open(str(data_file), 'r') as file:
        content = file.read()
    return content


def fetch_annotation(reference_id, reference_type=None):
    return _get_content('data/' + reference_id + '.gff3'), 'gff', 'ncbi'


def fetch_sequence(reference_id, reference_source=None):
    return _get_content('data/' + reference_id + '.sequence')


@pytest.mark.parametrize(
    'hgvs_description, normalized_description',
    [('NG_012337.1:g.4delins7_31',
      'NG_012337.1:g.4delins7_31'),
     ('NG_012337.1:g.26_31del',
      'NG_012337.1:g.29_34del'),
     ('NG_012337.1:g.4delins7_100',
      'NG_012337.1:g.[3_4insGGTT;5_6ins12_100]'),
     ('NG_012337.1:g.4delins7_50',
      'NG_012337.1:g.[3_4insGGTT;5_6ins12_50]'),
     ('NG_012337.1:g.26_31del',
      'NG_012337.1:g.29_34del'),
     ('NG_012337.1:g.4A>T',
      'NG_012337.1:g.4C>T'),
     ('NG_012337.1:g.100_200delins100_101',
      'NG_012337.1:g.102_200del'),
     # To be curated.
     # --------------
     # - Fix the ambiguity for the length/location.
     # ('NG_012337.1:g.100_200>400',
     #  'NG_012337.1:g.100_200delins100'),
     ])
def test_mutalyzer3(hgvs_description, normalized_description, monkeypatch):
    monkeypatch.setattr('retriever.retriever.fetch_annotations',
                        fetch_annotation)
    monkeypatch.setattr('retriever.retriever.fetch_sequence',
                        fetch_sequence)
    assert mutalyzer3(hgvs_description) == normalized_description
