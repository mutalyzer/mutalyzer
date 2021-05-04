import pytest

from normalizer.name_checker import name_check
from normalizer.reference import retrieve_reference, slice_to_selector
from normalizer.converter.to_rna import to_rna_reference_model, to_rna_variants
from normalizer.reference import get_selector_model
from mutalyzer_hgvs_parser import to_model

from .commons import code_in, patch_retriever

TESTS = [
    {
        "keywords": ["rna", "equal", "genomic"],
        "input": "NG_012337.1(NM_003002.2):r.275=",
        "normalized": "NG_012337.1(NM_003002.2):r.=",
        "to_test": True,
    },
    {
        "keywords": ["rna", "uncertain", "genomic"],
        "input": "NG_012337.1(NM_003002.2):r.?",
        "normalized": "NG_012337.1(NM_003002.2):r.?",
        "to_test": False,
    },
    {
        "keywords": ["rna", "substitution", "genomic"],
        "input": "NG_012337.1(NM_003002.2):r.275a>c",
        "normalized": "NG_012337.1(NM_003002.2):r.275a>c",
        "to_test": True,
    },
    {
        "keywords": ["rna", "substitution", "mRNA"],
        "input": "NM_003002.2:r.275a>c",
        "normalized": "NM_003002.2:r.275a>c",
        "to_test": True,
    },
    {
        "keywords": ["rna", "substitution", "ncRNA"],
        "input": "NR_038420.1:r.75g>a",
        "normalized": "NR_038420.1:r.75g>a",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "genomic"],
        "input": "NG_012337.1(NM_003002.2):r.273del",
        "normalized": "NG_012337.1(NM_003002.2):r.274del",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "mRNA"],
        "input": "NM_003002.2:r.273del",
        "normalized": "NM_003002.2:r.274del",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "ncRNA"],
        "input": "NR_038420.1:r.74del",
        "normalized": "NR_038420.1:r.75del",
        "to_test": True,
    },
]


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(t_type):
            output.append((test["input"], test[t_type]))
    return output


@pytest.mark.parametrize(
    "input_description, normalized",
    get_tests(TESTS, "normalized"),
)
def test_rna(input_description, normalized):
    d = name_check(input_description)
    assert d["normalized_description"] == normalized


@pytest.mark.parametrize(
    "r_id, s_id, expected",
    [
        (
            "TEST_REF",
            "NM_PLUS",
            {
                "annotations": {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 0},
                        "end": {"type": "point", "position": 636},
                    },
                    "id": "TEST_REF",
                    "type": "record",
                    "qualifiers": {
                        "name": "TEST_REF",
                        "chromosome": "TEST_CHR",
                        "mol_type": "genomic DNA",
                    },
                    "features": [
                        {
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": 0},
                                "end": {"type": "point", "position": 636},
                                "strand": 1,
                            },
                            "id": "PLUS",
                            "type": "gene",
                            "qualifiers": {"name": "PLUS"},
                            "features": [
                                {
                                    "id": "NM_PLUS",
                                    "type": "mRNA",
                                    "location": {
                                        "type": "range",
                                        "start": {"type": "point", "position": 0},
                                        "end": {"type": "point", "position": 636},
                                        "strand": 1,
                                    },
                                    "features": [
                                        {
                                            "id": "exon-NM_PLUS-1",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 54,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                        {
                                            "id": "exon-NM_PLUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 54,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 636,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                        {
                                            "id": "NP_PLUS",
                                            "type": "CDS",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 563,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                    ],
                                }
                            ],
                        }
                    ],
                },
                "sequence": {
                    "seq": "ATGGGGTCACAGTGGTTATAAAGAGTTATACCCTGAAGAATTTGAAACAGACAGTAGTGATCAGCAAGATATTACCAACGGGAAGAAAACATCTCCCCAGGTAAAGTCATCTACCCATGAATCCCGCAAACACAAGAAGTCAAAGAAATCCCACAAAAAAAAGCAGAAAAAAAGGTCACACAAAAAACAGAAGAAAAGCAAAAAGGAAGCCACAGATATAACAGCAGATTCCTCGAGTGAGTTCTCAGAAGAAACTGGGGCTTCTGGTACAAGGAAAGGGAAACAACCACATAAACGCAAGAAAAAATCCAGGAAAAAGTCTCTCAAAAAACCTGCTTTATTCTTAGAGGCAGAAAGTAACACTTCACATTCAGATGATTCAGCATCCAGCAGTTCTGAGGAAAGTGAGGAAAGAGACACTAAGAAAACCAAAAGGAAAAAGAGAGAGAAAAAAGCCCATACCTCTGTAGCCAACAATGAAATACAGGAGAGGACAAACAAACGCACAAATTGGAAAGTAGCTACAGATGAAAGGTCTGCTGAGAGCTCAGAGGATGACTAAATGGGAAACACTTTTGTTTTCCACATGACTGTGGATATTTACAGTTCTTACTCCTTGTGGTTTTGCCAGTGACT"
                },
            },
        ),
        (
            "TEST_REF",
            "NR_PLUS",
            {
                "annotations": {
                    "qualifiers": {
                        "name": "TEST_REF",
                        "chromosome": "TEST_CHR",
                        "mol_type": "genomic DNA",
                    },
                    "id": "TEST_REF",
                    "type": "record",
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 0},
                        "end": {"type": "point", "position": 503},
                    },
                    "features": [
                        {
                            "qualifiers": {"name": "PLUS"},
                            "id": "PLUS",
                            "type": "gene",
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": 0},
                                "end": {"type": "point", "position": 503},
                                "strand": 1,
                            },
                            "features": [
                                {
                                    "id": "NR_PLUS",
                                    "type": "ncRNA",
                                    "location": {
                                        "type": "range",
                                        "start": {"type": "point", "position": 0},
                                        "end": {"type": "point", "position": 503},
                                        "strand": 1,
                                    },
                                    "features": [
                                        {
                                            "id": "exon-NR_PLUS-1",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 101,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                        {
                                            "id": "exon-NR_PLUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 101,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 202,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                        {
                                            "id": "exon-NR_PLUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 202,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 503,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                    ],
                                }
                            ],
                        }
                    ],
                },
                "sequence": {
                    "seq": "TTGCTTTAACAATACATGTGATGTGTCATATTACAGATGGGGTCACAGTGGTTATAAAGAGTTATACCCTGAAGAATTTGAAACAGACAGGTAAGGAAAATGAAACAATAGGTGATTTAGTTAGGAGGTGATCCCTGTTTTCCCATTCTCTGGTTGATGTTTGGCATGTCTGTAAGCATTTTGGTTTTTATATATAGTATTCTTTTAACTTTTCTTATTAGTAGTGATCAGCAAGATATTACCAACGGGAAGAAAACATCTCCCCAGGTAAAGTCATCTACCCATGAATCCCGCAAACACAAGAAGTCAAAGAAATCCCACAAAAAAAAGCAGAAAAAAAGGTCACACAAAAAACAGAAGAAAAGCAAAAAGGAAGCCACAGATATAACAGCAGATTCCTCGAGTGAGTTCTCAGAAGAAACTGGGGCTTCTGGTACAAGGAAAGGGAAACAACCACATAAACGCAAGAAAAAATCCAGGAAAAAGTCTCTCAAAAAACCTGC"
                },
            },
        ),
        (
            "TEST_REF",
            "NM_MINUS",
            {
                "annotations": {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 0},
                        "end": {"type": "point", "position": 954},
                    },
                    "qualifiers": {
                        "name": "TEST_REF",
                        "chromosome": "TEST_CHR",
                        "mol_type": "genomic DNA",
                    },
                    "type": "record",
                    "id": "TEST_REF",
                    "features": [
                        {
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": 0},
                                "end": {"type": "point", "position": 954},
                                "strand": -1,
                            },
                            "qualifiers": {"name": "MINUS"},
                            "type": "gene",
                            "id": "MINUS",
                            "features": [
                                {
                                    "id": "NM_MINUS",
                                    "type": "mRNA",
                                    "location": {
                                        "type": "range",
                                        "start": {"type": "point", "position": 0},
                                        "end": {"type": "point", "position": 954},
                                        "strand": -1,
                                    },
                                    "features": [
                                        {
                                            "id": "exon-NM_MINUS-1",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 603,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 954,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NM_MINUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 502,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 603,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NM_MINUS-3",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 201,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 502,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NM_MINUS-4",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 201,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "NP_MINUS",
                                            "type": "CDS",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 125,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 855,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                    ],
                                }
                            ],
                        }
                    ],
                },
                "sequence": {
                    "seq": "TGTACTCCAAAGTCTTTCTAATGTTGCTTTAATTTCCAAAAATGTATGCATTGCTTTAACAATACATGTGATGTGTCATATTACAGATGGGGTCACAGTGGTTATAAAGAGTTATACCCTGAAGAATTTGAAACAGACAGGTAAGGAAAATAGGCTTACTGAAAGAAACTAAGATGGTACAAAATCTGTATTATAAATTGAGGTTGATGTTTGGCATGTCTGTAAGCATTTTGGTTTTTATATATAGTATTCCATACAGTAACTCAGTATGGCAGCTTAGAATTTTTACCTTCATTTTAAAGATGAGGAAACAAAAACTCAATGAGAATATTAAAGTGTTAAAGTATACATTAAAGTGCTTATTTAAAATTCAGATGTTAACCTCAATTTTTTAATCTAGAATGCAAAATATTAAAATAATACGCTTTTTTTTTACATAAAAGCTTCTATTTTTTAACTTTTCTTATTAGTAGTGATCAGCAAGATATTACCAACGGGAAGAAGTGAGTTCTCAGAAGAAACTGGGGCTTCTGGTACAAGGAAAGGGAAACAACCACATAAACGCAAGAAAAAATCCAGGAAAAAGTCTCTCAAAAAACCTGCGAAAAAGAGAGAGAAAAAAGCCCATACCTCTGTAGCCAACAATGAAATACAGGAGAGGACAAACAAACGCACAAATTGGAAAGTAGCTACAGATGAAAGGTCTGCTGAGAGCTCAGAGGATGACTAAATGGGAAACACTTTTGTTTTCCACATGACTGTGGATATTTACAGTTCTTACTCCTTGTGGTTTTGCCAGTGACTCTTGTTCAGCACGGGGCCTGAGGTCAGAGCTGTCTTGTGCCATCTGTCATTTCTGACAGACGTCTTGTCTTCTATTTTGGCGTTAAGCTTGATCCCCTTTTCTTGTTAAAAGGGAATCTGGTATTTTGTTATGAAGGTTTCTTGAAGAGA"
                },
            },
        ),
        (
            "TEST_REF",
            "NR_MINUS",
            {
                "annotations": {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 0},
                        "end": {"type": "point", "position": 993},
                    },
                    "id": "TEST_REF",
                    "qualifiers": {
                        "name": "TEST_REF",
                        "chromosome": "TEST_CHR",
                        "mol_type": "genomic DNA",
                    },
                    "type": "record",
                    "features": [
                        {
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": 0},
                                "end": {"type": "point", "position": 993},
                                "strand": -1,
                            },
                            "id": "MINUS",
                            "qualifiers": {"name": "MINUS"},
                            "type": "gene",
                            "features": [
                                {
                                    "id": "NR_MINUS",
                                    "type": "ncRNA",
                                    "location": {
                                        "type": "range",
                                        "start": {"type": "point", "position": 0},
                                        "end": {"type": "point", "position": 993},
                                        "strand": -1,
                                    },
                                    "features": [
                                        {
                                            "id": "exon-NR_MINUS-1",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 542,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 993,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NR_MINUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 191,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 542,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NR_MINUS-3",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 191,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                    ],
                                }
                            ],
                        }
                    ],
                },
                "sequence": {
                    "seq": "AGTCTTTCTAATGTTGCTTTAATTTCCAAAAATGTATGCATTGCTTTAACAATACATGTGATGTGTCATATTACAGATGGGGTCACAGTGGTTATAAAGAGTTATACCCTGAAGAATTTGAAACAGACAGGTAAGGAAAATAGGCTTACTGAAAGAAACTAAGATGGTACAAAATCTGTATTATAAATTGAGGTTGATGTTTGGCATGTCTGTAAGCATTTTGGTTTTTATATATAGTATTCCATACAGTAACTCAGTATGGCAGCTTAGAATTTTTACCTTCATTTTAAAGATGAGGAAACAAAAACTCAATGAGAATATTAAAGTGTTAAAGTATACATTAAAGTGCTTATTTAAAATTCAGATGTTAACCTCAATTTTTTAATCTAGAATGCAAAATATTAAAATAATACGCTTTTTTTTTACATAAAAGCTTCTATTTTTTAACTTTTCTTATTAGTAGTGATCAGCAAGATATTACCAACGGGAAGAAAACATCTCCCCAGGTAAAGTCATCTACCCATGAATCCCGCAAACACAAGCTTTATTCTTAGAGGCAGAAAGTAACACTTCACATTCAGATGATTCAGCATCCAGCAGTTCTGAGGAAAGTGAGGAAAGAGACACTAAGAAAACCAAAAGGAAAAAGAGAGAGAAAAAAGCCCATACCTCTGTAGCCAACAATGAAATACAGGAGAGGACAAACAAACGCACAAATTGGAAAGTAGCTACAGATGAAAGGTCTGCTGAGAGCTCAGAGGATGACTAAATGGGAAACACTTTTGTTTTCCACATGACTGTGGATATTTACAGTTCTTACTCCTTGTGGTTTTGCCAGTGACTCTTGTTCAGCACGGGGCCTGAGGTCAGAGCTGTCTTGTGCCATCTGTCATTTCTGACAGACGTCTTGTCTTCTATTTTGGCGTTAAGCTTGATCCCCTTTTCTTGTTAAAAGGGAATCTGGTATTTTGTTATGAAGGTTTCTTGAAGAGA"
                },
            },
        ),
    ],
)
def test_to_rna_reference_model(r_id, s_id, expected):

    model = retrieve_reference(r_id)
    rna_model = to_rna_reference_model(model, s_id)
    # import json
    # print(json.dumps(rna_model, indent=2))
    # print(len(rna_model["sequence"]["seq"]))
    assert rna_model == expected


def test_to_rna_variants():
    variants = to_model("[10_20delinsT;120_141delinsA;150_180delinsA;10_200delinsT;10_189delinsT;10_1200delinsT]", "variants")
    print(variants)
    sequences = {"TEST_REF": retrieve_reference("TEST_REF"),
                 "reference": retrieve_reference("TEST_REF")}
    s_model = get_selector_model(retrieve_reference("TEST_REF")["annotations"], "NM_PLUS")
    to_rna_variants(variants, sequences, s_model)
    assert 1 == 0
