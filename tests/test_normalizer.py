import pytest
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


@pytest.mark.parametrize(
    'hgvs_description, normalized_description',
    [('NG_012337.1:g.4delins7_31',
      'NG_012337.1:g.4delins7_31'),
     ('NG_012337.1:g.26_31del',
      'NG_012337.1:g.29_34del'),
     ('NG_017013.2:g.1_32772del',
      'NG_017013.2:g.1_32772del'),
     ('NG_017013.2:g.[16985A>T;17013_17014del]',
      'NG_017013.2:g.[16985A>T;17013_17014del]'),
     ('NG_017013.2:g.[16985A>T;17011_17012del]',
      'NG_017013.2:g.[16985A>T;17013_17014del]'),
     ('NG_017013.2:g.17013_17014del',
      'NG_017013.2:g.17013_17014del'),
     ('NG_017013.2:g.17011_17012del',
      'NG_017013.2:g.17013_17014del'),
     ('NG_017013.2:g.17415_17417delinsGCG',
      'NG_017013.2:g.[17415C>G;17417A>G]'),
     ('NG_017013.2:g.17496_17497insAGCTGCTCAGATAGCGA',
      'NG_017013.2:g.17496_17497ins[A;17481_17496]'),
     ('NG_012337.1:g.4delins7_50',
      'NG_012337.1:g.[3_4insGGTT;5_6ins12_50]'),
     ('NG_012337.1:g.26_31del',
      'NG_012337.1:g.29_34del'),
     ('NG_012337.1:g.4C>T',
      'NG_012337.1:g.4C>T'),
     ('NG_012337.1:g.100_200delins100_101',
      'NG_012337.1:g.102_200del'),
     ('NG_017013.2:g.19258dup',
      'NG_017013.2:g.19258dup'),
     ('NG_017013.2:g.16508_16509dup',
      'NG_017013.2:g.16508_16509dup'),
     ('NG_012337.1:g.90_91insC',
      'NG_012337.1:g.90_91insC'),
     ('NG_017013.2:g.17471_17471del',
      'NG_017013.2:g.17471del'),
     ('NG_017013.2:g.18748_18750delinsCAT',
      'NG_017013.2:g.18749G>A'),
     ('NG_017013.2:g.17394_17395insC',
      'NG_017013.2:g.17394dup'),
     ('NG_012337.1(NM_003002.2):c.274G>T',
      'NG_012337.1:g.7125G>T'),
     ('NG_012337.1(NM_003002.2):c.1del',
      'NG_012337.1:g.5062del'),
     ('NG_012337.1(NM_003002.2):c.-1del',
      'NG_012337.1:g.5061del'),
     ('NG_012337.1(NM_003002.2):c.52+1del',
      'NG_012337.1:g.5114del'),
     ('NG_012337.1(NM_003002.2):c.*824del',
      'NG_012337.1:g.13948del'),
     ('NG_012337.1(NM_003002.2):c.*824+10del',
      'NG_012337.1:g.13958del'),
     ('NM_003002.4:c.1del',
      'NM_003002.4:n.36del'),
     ('NG_029724.1(NM_004321.7):c.101del',
      'NG_029724.1:g.27557del'),
     ('NG_029724.1(NM_004321.7):c.101del',
      'NG_029724.1:g.27557del'),
     ('NG_029724.1:g.10_20del',
      'NG_029724.1:g.10_20del'),
     ('NG_029724.1:g.10_20del11',
      'NG_029724.1:g.10_20del'),
     ('NG_029724.1:g.10delG',
      'NG_029724.1:g.10del'),
     ('NG_029724.1:g.10_11delGT',
      'NG_029724.1:g.10_11del'),
     ('NG_008835.1(CDH23_v001):c.1449+846delA',
      'NG_008835.1:g.255529del'),
     ('NG_009113.2(NR2E3_v001):c.948delC',
      'NG_009113.2:g.8038del'),
     ('NG_007107.2(MECP2_v001):c.378-17delT',
      'NG_007107.2:g.110661del'),
     ('NG_009113.2(NR2E3_v002):c.948delC',
      'NG_009113.2:g.8038del'),
     ('NG_009113.2(NM_016346.4):c.948delC',
      'NG_009113.2:g.8038del'),
     ('NG_009113.2(NM_014249.4):c.948delC',
      'NG_009113.2:g.8038del'),
     ('NG_009497.1(NM_206933.2):c.8682-19dupT',
      'NG_009497.1:g.561208dup'),
     ('NG_009497.1(USH2A_v001):c.8682-19dup',
      'NG_009497.1:g.561208dup'),
     ])
def test_mutalyzer3(hgvs_description, normalized_description, monkeypatch):
    monkeypatch.setattr('retriever.retriever.fetch_annotations',
                        fetch_annotation)
    monkeypatch.setattr('retriever.retriever.fetch_sequence',
                        fetch_sequence)
    assert mutalyzer3(hgvs_description) == normalized_description


@pytest.mark.parametrize(
    'hgvs_description, exception_text',
    [
     ('NG_012337.1:c.10del',
      'No selector mentioned for c. with a genomic reference.'),
     ('NG_029724.1(NM_004321.2):c.101del',
      'No exons.'),
     ('NM_001003806.1:c.323delG',
      'No CDS.'),
     ('NM_003002.4:c.205-40dupC',
      'Intronic positions for mRNA reference.'),
     ('NG_029724.1:g.10delA',
      'Description sequence mismatch.'),
     ('NG_029724.1:g.10delGT',
      'Description sequence mismatch.'),
     ('NM_003002.4:c.10del5',
      'Description sequence mismatch.'),
     ('NM_003002.4:c.10_12insA',
      'Not consecutive positions in insertion.'),
     ('NM_003002.4:c.10_11ins6',
      'Length in inserted not supported.'),
     ('NM_003002.4:c.10_11ins(5)',
      'Length in inserted not supported.'),
     ('NM_003002.4:c.10_11ins[5]',
      'Length in inserted not supported.'),
     ('NM_003002.4:c.10_11ins[50_70;5]',
      'Length in inserted not supported.'),
     ('NG_012337.1:g.100_200>400',
      'Length in inserted not supported.'),
     ('NM_003002.4:c.10AAA[50]',
      'Variant type not supported.'),
     ('NG_029724.1(NM_004321.7):c.100_*900000del',
      'Out of range.'),
     ('NG_009497.1(KCTD3_v001):c.8682-19dup',
      'No exons.'),
     ('NG_029724.1:g.10_5del',
      'End position is smaller then start position.'),
     ('NM_003002.4:c.*1_10del',
      'End position is smaller then start position.'),
     ('NG_012337.1:g.10delins20_15',
      'End position is smaller then start position.'),
     ('NG_012337.1:g.10del20_15',
      'End position is smaller then start position.'),
     ])
def test_mutalyzer3_exceptions(hgvs_description, exception_text, monkeypatch):
    monkeypatch.setattr('retriever.retriever.fetch_annotations',
                        fetch_annotation)
    monkeypatch.setattr('retriever.retriever.fetch_sequence',
                        fetch_sequence)
    with pytest.raises(Exception) as exc:
        mutalyzer3(hgvs_description)
    assert exception_text == str(exc.value)
