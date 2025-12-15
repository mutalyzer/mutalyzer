import pytest

from mutalyzer.normalizer import normalize

from .commons import code_in, monkey_patches
from .variants_set import TESTS_ALL


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(t_type):
            output.append((test["input"], test[t_type]))
    return output


def get_coding_protein_equivalent(tests):
    output = []
    for test in tests:
        if test.get("to_test") and test.get("coding_protein_descriptions"):
            for pair in test["coding_protein_descriptions"]:
                if len(pair) == 2:
                    output.append(pair)
    return output


def get_tests_rna_protein(tests):
    output = []
    for test in tests:
        if test.get("rna_description") and test.get("protein_description"):
            output.append((test["rna_description"], test["protein_description"]))
    return output


@pytest.mark.parametrize(
    "input_description, normalized", get_tests(TESTS_ALL, "normalized")
)
def test_normalize(input_description, normalized):
    d = normalize(input_description)
    assert d["normalized_description"] == normalized


@pytest.mark.parametrize("input_description, genomic", get_tests(TESTS_ALL, "genomic"))
def test_genomic(input_description, genomic):
    d = normalize(input_description)
    if d["equivalent_descriptions"].get("g"):
        assert d["equivalent_descriptions"]["g"][0]["description"] == genomic


@pytest.mark.parametrize(
    "input_description, coding", get_tests(TESTS_ALL, "coding_protein_descriptions")
)
def test_coding(input_description, coding):
    d = normalize(input_description)
    coding = [c[0] for c in coding]
    if d["equivalent_descriptions"].get("c"):
        name_check_coding = [
            c["description"] for c in d["equivalent_descriptions"]["c"]
        ]
    assert set(coding).issubset(set(name_check_coding))


@pytest.mark.parametrize(
    "input_description, protein_description",
    get_tests(TESTS_ALL, "protein_description"),
)
def test_protein(input_description, protein_description):

    normalized_output = normalize(input_description)
    normalizer_protein = normalized_output["protein"]["description"]

    assert normalizer_protein == protein_description


@pytest.mark.parametrize(
    "coding_description, protein_description",
    get_coding_protein_equivalent(TESTS_ALL),
)
def test_protein_equivalent(coding_description, protein_description):
    normalized_output = normalize(coding_description)
    normalizer_protein = normalized_output["protein"]["description"]
    assert normalizer_protein == protein_description


@pytest.mark.parametrize(
    "rna_description, protein_description",
    get_tests_rna_protein(TESTS_ALL),
)
def test_rna_protein(rna_description, protein_description):

    normalized_output = normalize(rna_description)
    normalizer_protein = normalized_output["protein"]["description"]

    assert normalizer_protein == protein_description


@pytest.mark.parametrize(
    "input_description, rna_description",
    get_tests(TESTS_ALL, "rna_description"),
)
def test_rna(input_description, rna_description):

    normalized_output = normalize(input_description)
    normalized_rna = normalized_output["rna"]["description"]

    assert normalized_rna == rna_description


@pytest.mark.parametrize("input_description, codes", get_tests(TESTS_ALL, "errors"))
def test_errors(input_description, codes):
    assert [error["code"] for error in normalize(input_description)["errors"]] == codes


@pytest.mark.parametrize("input_description, codes", get_tests(TESTS_ALL, "infos"))
def test_infos(input_description, codes):
    assert [info["code"] for info in normalize(input_description)["infos"]] == codes


@pytest.mark.parametrize(
    "description, sequence, normalized",
    [
        ("1A>T", "A", "1A>T"),
        ("1A>T", "AA", "1A>T"),
        ("2del", "AAAT", "1_3A[2]"),
        ("[2del]", "AAAT", "1_3A[2]"),
        ("[1del;2del]", "AAAT", "1_3A[1]"),
        ("1_2insNG_012337.1:g.100", "AAAT", "1_2insT"),
        ("[9dup;14_15insCCTCT]", "CTCTCTCTCTCTCTTG", "10delins[CT[3];C]"),
        ("10delinsCTCTCTC", "CTCTCTCTCTCTCTTG", "10delins[CT[3];C]"),
        ("[1_2del;5_10inv]", "CACACCCCCA", "3_10delinsTGGGGG"),
        ("[2_7del;8_34inv;35_36del]", "GTTCGCGGGGAAAGGAAAAAAGCCGCCGGGCAGGAAA", "[2_22delinsCCTGCCCGGCGGCTTTTTT;25delinsT[3];28_30del;33_37del]"),
    ],
)
def test_only_variants(description, sequence, normalized):
    assert (
        normalize(description, True, sequence)["normalized_description"] == normalized
    )


@pytest.mark.parametrize(
    "description, sequence, codes",
    [
        ("1C>T", "A", ["ESEQUENCEMISMATCH"]),
        ("1C>T", "AA", ["ESEQUENCEMISMATCH"]),
        ("2C>T", "AA", ["ESEQUENCEMISMATCH"]),
        ("1delC", "A", ["ESEQUENCEMISMATCH"]),
        ("1delC", "AA", ["ESEQUENCEMISMATCH"]),
        ("2delC", "AA", ["ESEQUENCEMISMATCH"]),
    ],
)
def test_only_variants_errors(description, sequence, codes):
    assert codes == [
        error["code"] for error in normalize(description, True, sequence)["errors"]
    ]


@pytest.mark.parametrize(
    "description, errors",
    [
        (
            "NM_003002.4:c.50+10del",
            [
                {
                    "code": "EINTRONIC",
                    "positions": {
                        "NM_003002.4": [
                            "50+10"
                        ]
                    },
                    "details": "Intronic position `50+10` was used with a non intronic reference sequence `NM_003002.4`.",
                    "options": [
                        {
                            "assembly_id": "GRCh38",
                            "description": "NC_000011.10(NM_003002.4):c.50+10del"
                        },
                        {
                            "assembly_id": "GRCh37",
                            "description": "NC_000011.9(NM_003002.4):c.50+10del"
                        }
                    ],
                    "assemblies": {
                        "GRCh38": {
                            "NM_003002.4": {
                                "chr_id": "NC_000011.10",
                                "slices_differ": False
                            }
                        },
                        "GRCh37": {
                            "NM_003002.4": {
                                "chr_id": "NC_000011.9",
                                "slices_differ": False
                            }
                        }
                    }
                }
            ],
        ),
        (
            "NM_003002.4:c.50+10_50+12del",
            [
                {
                    "code": "EINTRONIC",
                    "positions": {
                        "NM_003002.4": [
                            "50+10",
                            "50+12"
                        ]
                    },
                    "details": "Intronic positions `50+10`, `50+12` were used with a non intronic reference sequence `NM_003002.4`.",
                    "options": [
                        {
                            "assembly_id": "GRCh38",
                            "description": "NC_000011.10(NM_003002.4):c.50+10_50+12del"
                        },
                        {
                            "assembly_id": "GRCh37",
                            "description": "NC_000011.9(NM_003002.4):c.50+10_50+12del"
                        }
                    ],
                    "assemblies": {
                        "GRCh38": {
                            "NM_003002.4": {
                                "chr_id": "NC_000011.10",
                                "slices_differ": False
                            }
                        },
                        "GRCh37": {
                            "NM_003002.4": {
                                "chr_id": "NC_000011.9",
                                "slices_differ": False
                            }
                        }
                    }
                }
            ],
        ),
        (
            "NG_012337.3(NM_003002.4):c.274delinsNM_003002.4(SDHD):c.310+10",
            [
                {
                    "code": "EINTRONIC",
                    "positions": {
                        "NM_003002.4": [
                            "310+10"
                        ]
                    },
                    "details": "Intronic position `310+10` was used with a non intronic reference sequence `NM_003002.4`.",
                    "options": [
                        {
                            "assembly_id": "GRCh38",
                            "description": "NG_012337.3(NM_003002.4):c.274delinsNC_000011.10(NM_003002.4):c.310+10"
                        },
                        {
                            "assembly_id": "GRCh37",
                            "description": "NG_012337.3(NM_003002.4):c.274delinsNC_000011.9(NM_003002.4):c.310+10"
                        }
                    ],
                    "assemblies": {
                        "GRCh38": {
                            "NM_003002.4": {
                                "chr_id": "NC_000011.10",
                                "slices_differ": False
                            }
                        },
                        "GRCh37": {
                            "NM_003002.4": {
                                "chr_id": "NC_000011.9",
                                "slices_differ": False
                            }
                        }
                    }
                }
            ],
        ),
        (
            "NM_024426.4:c.[52+5_52+10del;100delinsNM_003002.4:52+2;169+1delins52+10_52+15]",
            [
                {
                    "code": "EINTRONIC",
                    "positions": {
                        "NM_024426.4": [
                            "52+5",
                            "52+10",
                            "169+1",
                            "52+10",
                            "52+15"
                        ],
                        "NM_003002.4": [
                            "52+2"
                        ]
                    },
                    "details": "Intronic positions `52+5`, `52+10`, `169+1`, `52+10`, `52+15` were used with a non "
                               "intronic reference sequence `NM_024426.4`. Intronic position `52+2` was used with a "
                               "non intronic reference sequence `NM_003002.4`.",
                    "options": [
                        {
                            "assembly_id": "GRCh38",
                            "description": "NM_024426.4:c.[52+5_52+10del;100delinsNC_000011.10(NM_003002.4):c.52+2;169+1delins52+10_52+15]"
                        },
                        {
                            "assembly_id": "GRCh37",
                            "description": "NC_000011.9(NM_024426.4):c.[52+5_52+10del;100delinsNC_000011.9(NM_003002.4):c.52+2;169+1delins52+10_52+15]"
                        }
                    ],
                    "assemblies": {
                        "GRCh38": {
                            "NM_003002.4": {
                                "chr_id": "NC_000011.10",
                                "slices_differ": False
                            }
                        },
                        "GRCh37": {
                            "NM_024426.4": {
                                "chr_id": "NC_000011.9",
                                "slices_differ": False
                            },
                            "NM_003002.4": {
                                "chr_id": "NC_000011.9",
                                "slices_differ": False
                            }
                        }
                    }
                }
            ]
                ,
        ),
        (
            "NG_012337.3(NM_003002.4):c.274delins[NM_003002.4(SDHD):c.310+10;NM_024426.4:c.100+50_100+55]",
            [
                {
                    "code": "EINTRONIC",
                    "positions": {
                        "NM_003002.4": [
                            "310+10"
                        ],
                        "NM_024426.4": [
                            "100+50",
                            "100+55"
                        ]
                    },
                    "details": "Intronic position `310+10` was used with a non intronic reference sequence "
                               "`NM_003002.4`. Intronic positions `100+50`, `100+55` were "
                               "used with a non intronic reference sequence `NM_024426.4`.",
                    "options": [
                        {
                            "assembly_id": "GRCh38",
                            "description": "NG_012337.3(NM_003002.4):c.274delins[NC_000011.10(NM_003002.4):c.310+10;NM_024426.4:c.100+50_100+55]"
                        },
                        {
                            "assembly_id": "GRCh37",
                            "description": "NG_012337.3(NM_003002.4):c.274delins[NC_000011.9(NM_003002.4):c.310+10;NC_000011.9(NM_024426.4):c.100+50_100+55]"
                        }
                    ],
                    "assemblies": {
                        "GRCh38": {
                            "NM_003002.4": {
                                "chr_id": "NC_000011.10",
                                "slices_differ": False
                            }
                        },
                        "GRCh37": {
                            "NM_024426.4": {
                                "chr_id": "NC_000011.9",
                                "slices_differ": False
                            },
                            "NM_003002.4": {
                                "chr_id": "NC_000011.9",
                                "slices_differ": False
                            }
                        }
                    }
                }
            ],
        ),
    ],
)
def test_intronic(monkeypatch, description, errors):
    def _get_chromosome_from_selector(assembly_id, ref_id):
        if assembly_id == "GRCh38" and ref_id in ["NM_003002.4", "NM_024426.4"]:
            return "NC_000011.10"
        if assembly_id == "GRCh37" and ref_id in ["NM_003002.4", "NM_024426.4"]:
            return "NC_000011.9"
        return None

    def _sequences_differ(r, s):
        if r == "NC_000011.10" and s in ["NM_003002.4", "NM_024426.4"]:
            return False
        elif r == "NC_000011.9" and s in ["NM_003002.4", "NM_024426.4"]:
            return False
        return True

    monkeypatch.setattr(
        "mutalyzer.description.get_chromosome_from_selector",
        _get_chromosome_from_selector,
    )
    monkeypatch.setattr(
        "mutalyzer.description._slices_differ",
        _sequences_differ,
    )
    assert (sorted(normalize(description)["errors"], key=lambda e: (e.get('code', ''), str(e.get('paths', [])))) ==
            sorted(errors, key=lambda e: (e.get('code', ''), str(e.get('paths', [])))))



@pytest.mark.parametrize(
    "description, chr_id, gene, gene_suggestions, errors",
    [
        (
                "SDHD:c.52+65del",
                "NG_012337.3",  # We use this instead of NC_000011.10
                "SDHD",
                {
                    'NG_012337.3': [
                        {"id": "NM_003002.3"},
                        {"id": "NM_001276506.2"}
                    ]
                },
                [{
                    "code": 'EGENEASREFERENCEID',
                    "details": "`SDHD` is an invalid reference identifier; it appears to be a gene name.",
                    "gene": 'SDHD',
                    "chr_ids": ["NG_012337.3"],
                    "options": {"NG_012337.3": [{"description": "NG_012337.3(NM_003002.3):c.52+65del",
                                                 "transcript_id": "NM_003002.3",
                                                 "chromosome_id": "NG_012337.3"
                                                 },
                                                {"description": "NG_012337.3(NM_001276506.2):c.52+65del",
                                                 "transcript_id": "NM_001276506.2",
                                                 "chromosome_id": "NG_012337.3"}]},
                    "paths": [("reference", "id")]}]
        ),
        (
                "SDHD:c.274delinsNM_003002.4:52",
                "NG_012337.3",
                "SDHD",
                {
                    "NG_012337.3": [
                        {"id": "NM_003002.3", "tag": "RefSeq Select"},
                        {"id": "NM_001276506.2"}
                    ]
                },
                [{
                    "code": "EGENEASREFERENCEID",
                    "details": "`SDHD` is an invalid reference identifier; it appears to be a gene name.",
                    "gene": "SDHD",
                    "chr_ids": ["NG_012337.3"],
                    "options": {"NG_012337.3": [{"description": "NG_012337.3(NM_003002.3):c.274delinsNM_003002.4:52",
                                                 "transcript_id": "NM_003002.3",
                                                 "chromosome_id": "NG_012337.3",
                                                 "tag": "RefSeq Select"},
                                                {"description": "NG_012337.3(NM_001276506.2):c.274delinsNM_003002.4:52",
                                                 "transcript_id": "NM_001276506.2",
                                                 "chromosome_id": "NG_012337.3"}]},
                    "paths": [("reference", "id")]}]

        ),
        (
            "UNKNOWNGENE:c.100A>G",
            None,
            "UNKNOWNGENE",
            None,
            [{
                "code": "ERETR",
                "details": "`UNKNOWNGENE` could not be retrieved; it may be an invalid reference identifier.",
                "paths": [[("reference", "id")]]
            }]
        ),
    ],
)
def test_gene_as_reference_id_error(monkeypatch, description, chr_id, gene, gene_suggestions, errors):
    def _get_gene_suggestions(gene_name):
        if gene_name == gene:
            return gene_suggestions
        return None

    def get_assembly_from_chr_id(gene_suggestions):
        return chr_id

    monkeypatch.setattr(
        "mutalyzer.description.get_gene_suggestions",
        _get_gene_suggestions,
    )
    monkeypatch.setattr(
        "mutalyzer.description.get_assembly_from_chr_id",
        get_assembly_from_chr_id,
    )
    assert normalize(description)["errors"] == errors
