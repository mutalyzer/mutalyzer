import pytest
from mutalyzer_hgvs_parser import to_model

from normalizer.converter.to_rna import (
    _get_location_type,
    _trim_to_exons,
    to_rna_reference_model,
    to_rna_variants,
)
from normalizer.name_checker import name_check
from normalizer.reference import get_selector_model, retrieve_reference

from .commons import code_in, patch_retriever

TESTS = [
    {
        "keywords": ["rna", "equal", "genomic"],
        "input": "NG_012337.1(NM_003002.2):r.275=",
        "normalized": "NG_012337.1(NM_003002.2):r.275=",
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
                    "seq": "auggggucacagugguuauaaagaguuauacccugaagaauuugaaacagacaguagugaucagcaagauauuaccaacgggaagaaaacaucuccccagguaaagucaucuacccaugaaucccgcaaacacaagaagucaaagaaaucccacaaaaaaaagcagaaaaaaaggucacacaaaaaacagaagaaaagcaaaaaggaagccacagauauaacagcagauuccucgagugaguucucagaagaaacuggggcuucugguacaaggaaagggaaacaaccacauaaacgcaagaaaaaauccaggaaaaagucucucaaaaaaccugcuuuauucuuagaggcagaaaguaacacuucacauucagaugauucagcauccagcaguucugaggaaagugaggaaagagacacuaagaaaaccaaaaggaaaaagagagagaaaaaagcccauaccucuguagccaacaaugaaauacaggagaggacaaacaaacgcacaaauuggaaaguagcuacagaugaaaggucugcugagagcucagaggaugacuaaaugggaaacacuuuuguuuuccacaugacuguggauauuuacaguucuuacuccuugugguuuugccagugacu"
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
                    "seq": "uugcuuuaacaauacaugugaugugucauauuacagauggggucacagugguuauaaagaguuauacccugaagaauuugaaacagacagguaaggaaaaugaaacaauaggugauuuaguuaggaggugaucccuguuuucccauucucugguugauguuuggcaugucuguaagcauuuugguuuuuauauauaguauucuuuuaacuuuucuuauuaguagugaucagcaagauauuaccaacgggaagaaaacaucuccccagguaaagucaucuacccaugaaucccgcaaacacaagaagucaaagaaaucccacaaaaaaaagcagaaaaaaaggucacacaaaaaacagaagaaaagcaaaaaggaagccacagauauaacagcagauuccucgagugaguucucagaagaaacuggggcuucugguacaaggaaagggaaacaaccacauaaacgcaagaaaaaauccaggaaaaagucucucaaaaaaccugc"
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
                    "seq": "uguacuccaaagucuuucuaauguugcuuuaauuuccaaaaauguaugcauugcuuuaacaauacaugugaugugucauauuacagauggggucacagugguuauaaagaguuauacccugaagaauuugaaacagacagguaaggaaaauaggcuuacugaaagaaacuaagaugguacaaaaucuguauuauaaauugagguugauguuuggcaugucuguaagcauuuugguuuuuauauauaguauuccauacaguaacucaguauggcagcuuagaauuuuuaccuucauuuuaaagaugaggaaacaaaaacucaaugagaauauuaaaguguuaaaguauacauuaaagugcuuauuuaaaauucagauguuaaccucaauuuuuuaaucuagaaugcaaaauauuaaaauaauacgcuuuuuuuuuacauaaaagcuucuauuuuuuaacuuuucuuauuaguagugaucagcaagauauuaccaacgggaagaagugaguucucagaagaaacuggggcuucugguacaaggaaagggaaacaaccacauaaacgcaagaaaaaauccaggaaaaagucucucaaaaaaccugcgaaaaagagagagaaaaaagcccauaccucuguagccaacaaugaaauacaggagaggacaaacaaacgcacaaauuggaaaguagcuacagaugaaaggucugcugagagcucagaggaugacuaaaugggaaacacuuuuguuuuccacaugacuguggauauuuacaguucuuacuccuugugguuuugccagugacucuuguucagcacggggccugaggucagagcugucuugugccaucugucauuucugacagacgucuugucuucuauuuuggcguuaagcuugauccccuuuucuuguuaaaagggaaucugguauuuuguuaugaagguuucuugaagaga"
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
                    "seq": "agucuuucuaauguugcuuuaauuuccaaaaauguaugcauugcuuuaacaauacaugugaugugucauauuacagauggggucacagugguuauaaagaguuauacccugaagaauuugaaacagacagguaaggaaaauaggcuuacugaaagaaacuaagaugguacaaaaucuguauuauaaauugagguugauguuuggcaugucuguaagcauuuugguuuuuauauauaguauuccauacaguaacucaguauggcagcuuagaauuuuuaccuucauuuuaaagaugaggaaacaaaaacucaaugagaauauuaaaguguuaaaguauacauuaaagugcuuauuuaaaauucagauguuaaccucaauuuuuuaaucuagaaugcaaaauauuaaaauaauacgcuuuuuuuuuacauaaaagcuucuauuuuuuaacuuuucuuauuaguagugaucagcaagauauuaccaacgggaagaaaacaucuccccagguaaagucaucuacccaugaaucccgcaaacacaagcuuuauucuuagaggcagaaaguaacacuucacauucagaugauucagcauccagcaguucugaggaaagugaggaaagagacacuaagaaaaccaaaaggaaaaagagagagaaaaaagcccauaccucuguagccaacaaugaaauacaggagaggacaaacaaacgcacaaauuggaaaguagcuacagaugaaaggucugcugagagcucagaggaugacuaaaugggaaacacuuuuguuuuccacaugacuguggauauuuacaguucuuacuccuugugguuuugccagugacucuuguucagcacggggccugaggucagagcugucuugugccaucugucauuucugacagacgucuugucuucuauuuuggcguuaagcuugauccccuuuucuuguuaaaagggaaucugguauuuuguuaugaagguuucuugaagaga"
                },
            },
        ),
    ],
)
def test_to_rna_reference_model(r_id, s_id, expected):
    model = retrieve_reference(r_id)
    rna_model = to_rna_reference_model(model, s_id)
    assert rna_model == expected


@pytest.mark.parametrize(
    "variants, expected",
    [
        (  # same intron
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 10},
                        "end": {"type": "point", "position": 20},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                    "inserted": [{"sequence": "T", "source": "description"}],
                },
            ],
            [],
        ),
        (  # intron exon with insertion
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 120},
                        "end": {"type": "point", "position": 141},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                    "inserted": [{"sequence": "A", "source": "description"}],
                },
            ],
            [],
        ),
        (  # intron exon without insertion
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 120},
                        "end": {"type": "point", "position": 141},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                },
            ],
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 135},
                        "end": {"type": "point", "position": 141},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                }
            ],
        ),
        (  # same exon
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 150},
                        "end": {"type": "point", "position": 180},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                    "inserted": [{"sequence": "A", "source": "description"}],
                },
            ],
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 150},
                        "end": {"type": "point", "position": 180},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                    "inserted": [{"sequence": "A", "source": "description"}],
                }
            ],
        ),
        (  # intron intron with insertion
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 10},
                        "end": {"type": "point", "position": 200},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                    "inserted": [{"sequence": "T", "source": "description"}],
                },
            ],
            [],
        ),
        (  # intron intron without insertion
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 10},
                        "end": {"type": "point", "position": 200},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                },
            ],
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 135},
                        "end": {"type": "point", "position": 189},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                }
            ],
        ),
        (
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 10},
                        "end": {"type": "point", "position": 189},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                    "inserted": [{"sequence": "T", "source": "description"}],
                },
            ],
            [],
        ),
        (  # intron intron over exon with insertion
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 10},
                        "end": {"type": "point", "position": 1200},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                    "inserted": [{"sequence": "T", "source": "description"}],
                },
            ],
            [],
        ),
        (  # intron intron over exon without insertion
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 10},
                        "end": {"type": "point", "position": 1201},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                },
            ],
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 135},
                        "end": {"type": "point", "position": 1200},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                }
            ],
        ),
        (  # exon exon over intron without insertion
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 150},
                        "end": {"type": "point", "position": 1100},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                },
            ],
            [
                {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 150},
                        "end": {"type": "point", "position": 1100},
                    },
                    "type": "deletion_insertion",
                    "source": "reference",
                }
            ],
        ),
    ],
)
def test_trim_to_exons(variants, expected):
    sequences = {
        "TEST_REF": retrieve_reference("TEST_REF"),
        "reference": retrieve_reference("TEST_REF"),
    }
    selector_model = get_selector_model(
        retrieve_reference("TEST_REF")["annotations"], "NM_PLUS"
    )
    exons = [e for exon in selector_model["exon"] for e in exon]
    # exons: [135, 189, 618, 1200]
    assert _trim_to_exons(variants, exons, sequences) == expected


@pytest.mark.parametrize(
    "location, exons, location_type",
    [
        ("134_134", [135, 189, 618, 1200], "same intron"),
        ("134_135", [135, 189, 618, 1200], "same intron"),
        ("134_136", [135, 189, 618, 1200], "adjacent intron exon"),
        ("134_188", [135, 189, 618, 1200], "adjacent intron exon"),
        ("134_189", [135, 189, 618, 1200], "adjacent intron exon"),
        ("134_190", [135, 189, 618, 1200], "intron intron"),
        ("134_617", [135, 189, 618, 1200], "intron intron"),
        ("134_618", [135, 189, 618, 1200], "intron intron"),
        ("134_619", [135, 189, 618, 1200], "intron exon"),
        ("134_1199", [135, 189, 618, 1200], "intron exon"),
        ("134_1200", [135, 189, 618, 1200], "intron exon"),
        ("134_1201", [135, 189, 618, 1200], "intron intron"),
        ("135_135", [135, 189, 618, 1200], "same exon"),
        ("135_136", [135, 189, 618, 1200], "same exon"),
        ("135_188", [135, 189, 618, 1200], "same exon"),
        ("135_189", [135, 189, 618, 1200], "same exon"),
        ("135_190", [135, 189, 618, 1200], "adjacent exon intron"),
        ("135_617", [135, 189, 618, 1200], "adjacent exon intron"),
        ("135_618", [135, 189, 618, 1200], "adjacent exon intron"),
        ("135_619", [135, 189, 618, 1200], "exon exon"),
        ("136_188", [135, 189, 618, 1200], "same exon"),
        ("136_189", [135, 189, 618, 1200], "same exon"),
        ("136_190", [135, 189, 618, 1200], "adjacent exon intron"),
        ("136_617", [135, 189, 618, 1200], "adjacent exon intron"),
        ("136_618", [135, 189, 618, 1200], "adjacent exon intron"),
        ("136_619", [135, 189, 618, 1200], "exon exon"),
        ("188_188", [135, 189, 618, 1200], "same exon"),
        ("188_189", [135, 189, 618, 1200], "same exon"),
        ("188_190", [135, 189, 618, 1200], "adjacent exon intron"),
        ("189_189", [135, 189, 618, 1200], "same intron"),
        ("617_617", [135, 189, 618, 1200], "same intron"),
        ("617_618", [135, 189, 618, 1200], "same intron"),
        ("617_619", [135, 189, 618, 1200], "adjacent intron exon"),
        ("618_618", [135, 189, 618, 1200], "same exon"),
        ("618_619", [135, 189, 618, 1200], "same exon"),
        ("618_1199", [135, 189, 618, 1200], "same exon"),
        ("618_1200", [135, 189, 618, 1200], "same exon"),
        ("618_1201", [135, 189, 618, 1200], "adjacent exon intron"),
        ("1199_1199", [135, 189, 618, 1200], "same exon"),
        ("1199_1200", [135, 189, 618, 1200], "same exon"),
        ("1199_1201", [135, 189, 618, 1200], "adjacent exon intron"),
        ("1200_1201", [135, 189, 618, 1200], "same intron"),
        ("188_188", [135, 189, 189, 1200], "same exon"),
        ("188_189", [135, 189, 189, 1200], "same exon"),
        ("189_189", [135, 189, 189, 1200], "same exon"),
        ("188_190", [135, 189, 189, 1200], "exon exon"),  # TODO: adjacent exon?
    ],
)
def test_get_location_type(location, exons, location_type):
    location = to_model(location, "location")
    assert _get_location_type(location, exons) == location_type
