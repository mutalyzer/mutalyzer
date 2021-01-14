import json
from pathlib import Path

import pytest

from normalizer.name_checker import name_check

from .test_set import TESTS_ALL

TEST_SET = [
    {
        "keywords": [
            "M2: deletion_in_frame",
            "Simple in-frame deletion should give a simple description on protein level.",
            "To be adapted: switch from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_058195.3",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_163del",
        "input": "NG_007485.1(NM_058195.3):c.?",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.?",
                # "AL449423.14(CDKN2A_i001):p.(Met54_Gly55delinsSer)",
                "NG_007485.1(NP_478102.2):p.?",
            ),
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: insertion_in_frame",
            "Simple in-frame insertion should give a simple description on protein level.",
            "To be adapted.",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162insATC",
        "input": "NG_007485.1(NM_058195.3):c.?",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.?",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsIleSer)"
                "NG_007485.1(NP_478102.2):p.?",
            )
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: deletion_insertion_in_frame",
            "Simple in-frame deletion/insertion should give a simple description on protein level.",
            "To be adapted.",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delinsATCCC",
        "input": "NG_007485.1(NM_058195.3):c.?",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.?",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_478102.2):p.?",
            )
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: deletion_insertion_list_in_frame",
            "Simple in-frame deletion-insertion of a list should give a simple description on protein level.",
            "To be adapted.",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delins[ATCCC]",
        "input": "NG_007485.1(NM_058195.3):c.?",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.?",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_478102.2):p.?",
            )
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: deletion_insertion_in_frame_complete",
            "Simple in-frame deletion-insertion of a list should give a simple description on protein level.",
            "To be adapted.",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delTGinsATCCC",
        "input": "NG_007485.1(NM_058195.3):c.?",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.?",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_478102.2):p.?",
            )
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: deletion_insertion_list_in_frame_complete",
            "Simple in-frame deletion-insertion of a list should give a simple "
            "description on protein level, also with the optional deleted "
            "sequence argument.",
            "To be adapted.",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delTGins[ATCCC]",
        "input": "NG_007485.1(NM_058195.3):c.?",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.?",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_478102.2):p.?",
            )
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: del_exon_unknown_offsets",
            "Deletion of an entire exon with unknown offsets should be possible.",
            "To be implemented.",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.632-?_681+?del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(BRCA2_v001):c.632-?_681+?del",
                "NG_012772.1(BRCA2_i001):p.(Val211Glufs*10)",
            )
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: delins_with_length",
            "Delins with explicit length of deleted sequence (bug #108).",
        ],
        "input": "NM_000193.2:c.108_109del2insG",
        "coding_protein_descriptions": {
            (
                "NM_000193.2(NM_000193.2):c.108_109delinsG",
                "NM_000193.2(NP_000184.1):p.(Lys38Serfs*2)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: fs_no_stop",
            "Frame shift yielding no stop codon should be described with "
            "uncertainty of the stop codon."
            "http://www.hgvs.org/mutnomen/FAQ.html#nostop",
        ],
        "input": "NM_001199.3:c.2188dup",
        "coding_protein_descriptions": {
            (
                "NM_001199.3(NM_001199.3):c.2188dup",
                "NM_001199.3(NP_001190.1):p.(Gln730Profs*?)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: ext_no_stop",
            "Extension yielding no stop codon should be described with "
            "uncertainty of the stop codon."
            "http://www.hgvs.org/mutnomen/FAQ.html#nostop",
            "To be implemented.",
        ],
        "input": "NM_000193.2:c.1388G>C",
        "coding_protein_descriptions": {
            (
                "NM_000193.2(NM_000193.2):c.1388G>C",
                "NM_000193.2(NP_000184.1):p.(*463Serext*?)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: ",
            "",
            "To be implemented.",
        ],
        "input": "",
        "coding_protein_descriptions": {
            (
                "",
                "",
            )
        },
        "to_test": False,
    },
]


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


def get_tests(tests):
    output = []
    for test in tests:
        if test.get("to_test") and test["coding_protein_descriptions"]:
            output.append((test["input"], test["coding_protein_descriptions"]))
    return output


@pytest.mark.parametrize(
    "input_description, coding_protein_descriptions", get_tests(TEST_SET)
)
def test_normalizer(input_description, coding_protein_descriptions, monkeypatch):
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", retrieve_raw)
    monkeypatch.setattr("normalizer.util.configuration", lambda: None)

    normalized_output = name_check(input_description)
    normalizer_descriptions = set(normalized_output["equivalent_descriptions"]["c"])

    assert coding_protein_descriptions.issubset(normalizer_descriptions)
