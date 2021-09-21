import pytest

from normalizer.algebra import _get_id, _get_hgvs, compare
from normalizer.reference import retrieve_reference

from .commons import code_in, patch_retriever


@pytest.mark.parametrize(
    "reference_id",
    ["NM_012459.2", "LRG_24"],
)
def test_get_id(reference_id):
    assert _get_id(reference_id) == {
        "input": reference_id,
        "type": "id",
        "reference": {"id": reference_id},
        "sequence": retrieve_reference(reference_id)["sequence"]["seq"],
    }


@pytest.mark.parametrize(
    "reference_id",
    ["NO_REF"],
)
def test_get_id_no_ref(reference_id):
    assert _get_id(reference_id) == {
        "input": reference_id,
        "type": "id",
        "reference": {"id": reference_id},
        "errors": [
            {
                "code": "ERETR",
                "details": "Reference NO_REF could not be retrieved.",
                "paths": [[]],
            }
        ],
    }


@pytest.mark.parametrize(
    "description, expected",
    [
        (
            "NM_012459.2:c.1del",
            {
                "input": "NM_012459.2:c.1del",
                "type": "hgvs",
                "reference": {"id": "NM_012459.2"},
                "sequence": {
                    "seq": "AAGTCGAGAGGCGGTGCACACCCGTCGCGCTGCGCAAACACAGCTGTCGGA"
                    "AGGTGGCGAGCCTGAGGCGAACAATGGCGGAGCTGGGCGAAGCCGATGAAGCGGAGTT"
                    "GCAGCGCCTGGTGGCCGCCGAGCAGCAGAAGGCGCAGTTTACTGCACAGGTGCATCAC"
                    "TTCATGGAGTTATGTTGGGATAAATGTGTGGAGAAGCCAGGGAATCGCCTAGACTCTC"
                    "GCACTGAAAATTGTCTCTCCAGCTGTGTAGACCGCTTCATTGACACCACTCTTGCCAT"
                    "CACCAGTCGGTTTGCCCAGATTGTACAGAAAGGAGGGCAGTAGGCCATCCCCCAGGAG"
                    "AATGACAGAAGCAAAGGACTTGTTACTAAGCAGATTTAAGGGTCAGTGGGGGAAGGCT"
                    "ATCAACCCATTGTCAGATCAGCATCAGGCTGTTATCAAGTCTGTTGGTGCTAAAAAGT"
                    "AAAAGATGAAATGTTCAAAGAGTGAAATTTATTTATTTGGAATTCAGAAATTCCAGGT"
                    "TGTATGACATCAGTTACTCAATAAGTGTGAATTCTCCAACTCTTCTTTTAATCCCATT"
                    "TTAGAATTTAATATAGAGATCTCTGATTGGCAGGAACACTAGAAATAAATGTTCCATG"
                    "GCCAGTAGTGCAAATGGGGGATTGTAGGTTTTGAAAAACCACCCTAAGCCATATTAAG"
                    "GGGGTTGGAAGAACCATCGAAGCCTAAGGCATAGAAGAAAATTTGGGGTTAAGAAAGA"
                    "TGAAGAACAAAAAACAGCTTTATTGCTTATACATGACCAAGAAAAGGAAAACATGGCA"
                    "AAAAAAAAAAAAAAAAA"
                },
            },
        )
    ],
)
def test_get_hgvs(description, expected):
    assert _get_hgvs(description) == expected


@pytest.mark.parametrize(
    "params, expected",
    [
        (
            ("AAAAA", "sequence", "ATAAAAA", "sequence", "2_3insT", "variant"),
            {"relation": "disjoint"},
        ),
        (
            (
                None,
                None,
                "LRG_24:g.5525_5532del",
                "hgvs",
                "LRG_24:g.5525_5533del",
                "hgvs",
            ),
            {"relation": "contains"},
        ),
        (
            (
                None,
                None,
                "LRG_24:g.5525_5532del",
                "hgvs",
                "LRG_24:g.5525_5533del",
                "hgvs",
            ),
            {"relation": "is_contained"},
        ),
        (
            (
                "LRG_24",
                "id",
                "LRG_24:g.5525_5532del",
                "hgvs",
                "LRG_24:g.5525_5533del",
                "hgvs",
            ),
            {"relation": "is_contained"},
        ),
        (
            (
                "LRG_24",
                "id",
                "LRG_24:g.5525_5532del",
                "hgvs",
                "LRG_24:g.5525_5533del",
                "hgvs",
            ),
            {"relation": "is_contained"},
        ),
    ],
)
def test_compare(params, expected):
    assert compare(*params) == expected
