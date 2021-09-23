import pytest

from normalizer.algebra import _get_hgvs_and_variant, _get_id, compare
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
                "sequence": "AAGTCGAGAGGCGGTGCACACCCGTCGCGCTGCGCAAACACAGCTGTCGG"
                "AAGGTGGCGAGCCTGAGGCGAACAATGGCGGAGCTGGGCGAAGCCGATGAAGCGGAGTTGCA"
                "GCGCCTGGTGGCCGCCGAGCAGCAGAAGGCGCAGTTTACTGCACAGGTGCATCACTTCATGG"
                "AGTTATGTTGGGATAAATGTGTGGAGAAGCCAGGGAATCGCCTAGACTCTCGCACTGAAAAT"
                "TGTCTCTCCAGCTGTGTAGACCGCTTCATTGACACCACTCTTGCCATCACCAGTCGGTTTGC"
                "CCAGATTGTACAGAAAGGAGGGCAGTAGGCCATCCCCCAGGAGAATGACAGAAGCAAAGGAC"
                "TTGTTACTAAGCAGATTTAAGGGTCAGTGGGGGAAGGCTATCAACCCATTGTCAGATCAGCA"
                "TCAGGCTGTTATCAAGTCTGTTGGTGCTAAAAAGTAAAAGATGAAATGTTCAAAGAGTGAAA"
                "TTTATTTATTTGGAATTCAGAAATTCCAGGTTGTATGACATCAGTTACTCAATAAGTGTGAA"
                "TTCTCCAACTCTTCTTTTAATCCCATTTTAGAATTTAATATAGAGATCTCTGATTGGCAGGA"
                "ACACTAGAAATAAATGTTCCATGGCCAGTAGTGCAAATGGGGGATTGTAGGTTTTGAAAAAC"
                "CACCCTAAGCCATATTAAGGGGGTTGGAAGAACCATCGAAGCCTAAGGCATAGAAGAAAATT"
                "TGGGGTTAAGAAAGATGAAGAACAAAAAACAGCTTTATTGCTTATACATGACCAAGAAAAGG"
                "AAAACATGGCAAAAAAAAAAAAAAAAAA",
                "reference_sequence": "AAGTCGAGAGGCGGTGCACACCCGTCGCGCATGCGCAAAC"
                "ACAGCTGTCGGAAGGTGGCGAGCCTGAGGCGAACAATGGCGGAGCTGGGCGAAGCCGATGA"
                "AGCGGAGTTGCAGCGCCTGGTGGCCGCCGAGCAGCAGAAGGCGCAGTTTACTGCACAGGTG"
                "CATCACTTCATGGAGTTATGTTGGGATAAATGTGTGGAGAAGCCAGGGAATCGCCTAGACT"
                "CTCGCACTGAAAATTGTCTCTCCAGCTGTGTAGACCGCTTCATTGACACCACTCTTGCCAT"
                "CACCAGTCGGTTTGCCCAGATTGTACAGAAAGGAGGGCAGTAGGCCATCCCCCAGGAGAAT"
                "GACAGAAGCAAAGGACTTGTTACTAAGCAGATTTAAGGGTCAGTGGGGGAAGGCTATCAAC"
                "CCATTGTCAGATCAGCATCAGGCTGTTATCAAGTCTGTTGGTGCTAAAAAGTAAAAGATGA"
                "AATGTTCAAAGAGTGAAATTTATTTATTTGGAATTCAGAAATTCCAGGTTGTATGACATCA"
                "GTTACTCAATAAGTGTGAATTCTCCAACTCTTCTTTTAATCCCATTTTAGAATTTAATATA"
                "GAGATCTCTGATTGGCAGGAACACTAGAAATAAATGTTCCATGGCCAGTAGTGCAAATGGG"
                "GGATTGTAGGTTTTGAAAAACCACCCTAAGCCATATTAAGGGGGTTGGAAGAACCATCGAA"
                "GCCTAAGGCATAGAAGAAAATTTGGGGTTAAGAAAGATGAAGAACAAAAAACAGCTTTATT"
                "GCTTATACATGACCAAGAAAAGGAAAACATGGCAAAAAAAAAAAAAAAAAA",
            },
        )
    ],
)
def test_get_hgvs(description, expected):
    assert _get_hgvs_and_variant(description) == expected


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
            {"relation": "is_contained"},
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
        (
            (
                None,
                None,
                "NG_008376.4:g.5525_5532del",
                "hgvs",
                "LRG_303:g.5525_5532del",
                "hgvs",
            ),
            {"relation": "equivalent"},
        ),
        (
            (
                None,
                None,
                "NG_012337.3:g.274>A",
                "hgvs",
                "274del",
                "variant",
            ),
            {"relation": "contains"},
        ),
        (
            (
                None,
                None,
                "NG_012337.3:g.274>A",
                "hgvs",
                "274del",
                "variant",
            ),
            {"relation": "contains"},
        ),
    ],
)
def test_compare(params, expected):
    assert compare(*params) == expected


@pytest.mark.parametrize(
    "params, expected",
    [
        (
            (
                None,
                "hgvs",
                "NG_012337.3:g.274>A",
                "HGVS",
                "1delete",
                "varianT",
            ),
            {
                "errors": {
                    "reference_type": [
                        {
                            "code": "EINVALIDINPUT",
                            "details": "'hgvs' not valid.",
                            "options": ["sequence", "id"],
                        }
                    ],
                    "lhs_type": [
                        {
                            "code": "EINVALIDINPUT",
                            "details": "'HGVS' not valid.",
                            "options": ["sequence", "variant", "hgvs"],
                        }
                    ],
                    "rhs_type": [
                        {
                            "code": "EINVALIDINPUT",
                            "details": "'varianT' not valid.",
                            "options": ["sequence", "variant", "hgvs"],
                        }
                    ],
                }
            },
        ),
        (
            (
                None,
                None,
                "NG_012337.3:g.274>A",
                "hgvs",
                "1delete",
                "variant",
            ),
            {
                "errors": {
                    "rhs": [
                        {
                            "code": "ESYNTAXUEOF",
                            "details": "Unexpected end of input.",
                            "pos_in_stream": 7,
                            "unexpected_character": "e",
                            "description": "1delete",
                            "expecting": [
                                "'(' for an uncertainty start or before a selector ID",
                                "':' between the reference part and the coordinate system",
                                "a reference / selector ID",
                            ],
                        }
                    ]
                }
            },
        ),
    ],
)
def test_compare_errors(params, expected):
    # import json
    #
    # print(json.dumps(compare(*params), indent=2))
    assert compare(*params) == expected
