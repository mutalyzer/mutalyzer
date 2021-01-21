import json
from pathlib import Path

import pytest

from normalizer.name_checker import name_check

TEST_SET = [
    {
        "keywords": [
            "M2: deletion_in_frame",
            "Simple in-frame deletion should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_163del",
        "input": "NG_007485.1(NM_000077.4):c.161_163del",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_163del",
                # "AL449423.14(CDKN2A_i001):p.(Met54_Gly55delinsSer)",
                "NG_007485.1(NP_000068.1):p.(Met54_Gly55delinsSer)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_206del",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGlu)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: insertion_in_frame",
            "Simple in-frame insertion should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162insATC",
        "input": "NG_007485.1(NM_000077.4):c.161_162insATC",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162insATC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsIleSer)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsIleSer)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205insATC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69insIle)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: deletion_insertion_in_frame",
            "Simple in-frame deletion/insertion should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delinsATCCC",
        "input": "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: deletion_insertion_list_in_frame",
            "Simple in-frame deletion-insertion of a list should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delins[ATCCC]",
        "input": "NG_007485.1(NM_000077.4):c.161_162delins[ATCCC]",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: deletion_insertion_in_frame_complete",
            "Simple in-frame deletion-insertion of a list should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delTGinsATCCC",
        "input": "NG_007485.1(NM_000077.4):c.161_162delTGinsATCCC",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: deletion_insertion_list_in_frame_complete",
            "Simple in-frame deletion-insertion of a list should give a simple "
            "description on protein level, also with the optional deleted "
            "sequence argument.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delTGins[ATCCC]",
        "input": "NG_007485.1(NM_000077.4):c.161_162delTGins[ATCCC]",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "to_test": True,
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
            "M2: fs_ext_no_stop",
            "Extension yielding no stop codon should be described with "
            "uncertainty of the stop codon."
            "http://www.hgvs.org/mutnomen/FAQ.html#nostop",
            "To be implemented.",
        ],
        "input": "NM_000193.2:c.1388_1389insC",
        "coding_protein_descriptions": {
            (
                "NM_000193.2(NM_000193.2):c.1388_1389insC",
                "NM_000193.2(NP_000184.1):p.(*463Cysext*?)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: synonymous_p_is",
            "Synonymous mutation should yield a p.(=) description.",
            "Switched from AB026906.1 to NG_012337.1 and from SDHD_v001 to NM_003002.2",
        ],
        "input": "NG_012337.1(NM_003002.2):c.276C>T",
        "coding_protein_descriptions": {
            (
                # "AB026906.1(SDHD_v001):c.276C>T",
                "NG_012337.1(NM_003002.2):c.276C>T",
                "NG_012337.1(NP_002993.1):p.(=)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: synonymous_p_is_alt_start",
            "Synonymous mutation should yield a p.(=) description, also with an "
            "alternative start codon.",
        ],
        "input": "NM_024426.4:c.1107A>G",
        "coding_protein_descriptions": {
            (
                "NM_024426.4(NM_024426.4):c.1107A>G",
                "NM_024426.4(NP_077744.3):p.(=)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: start_codon",
            "Mutation of start codon should yield a p.? description.",
            "Used NG_012337.1 instead of AB026906.1." "To be implemented.",
        ],
        # "input": "AB026906.1:c.1A>G",
        "input": "NG_012337.1(NM_003002.2):c.1A>G",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_003002.2):c.1A>G",
                "NG_012337.1(NP_002993.1):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: start_codon_alt_start",
            "Mutation of start codon should yield a p.? description, also with "
            "an alternative start codon.",
            "To be implemented.",
        ],
        "input": "NM_024426.4:c.1C>G",
        "coding_protein_descriptions": {
            (
                "NM_024426.4(NM_024426.4):c.1C>G",
                "NM_024426.4(NP_077744.3):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: start_codon_yield_start_p_is",
            "Silent mutation creating new start codon should yield a p.? "
            "description. The visualisation should also render the case for "
            "the new start codon.",
            "Used NG_012337.1 instead of AB026906.1." "To be implemented.",
        ],
        # "input": "AB026906.1:c.1A>T", # yields TTG start codon
        "input": "NG_012337.1(NM_003002.2):c.1A>T",  # yields TTG start codon
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_003002.2):c.1A>T",
                "NG_012337.1(NP_002993.1):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: start_codon_alt_start_yield_start_p_is",
            "Silent mutation creating new start codon should yield a p.? "
            "description, also with an alternative start codon. The "
            "visualisation should also render the case for the new start codon.",
            "To be implemented.",
        ],
        "input": "NM_024426.4:c.1C>A",  # yields ATG start codon
        "coding_protein_descriptions": {
            (
                "NM_024426.4(NM_024426.4):c.1C>A",
                "NM_024426.4(NP_077744.3):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: start_codon_yield_start",
            "Mutation creating new start codon should yield a p.? description. "
            "The visualisation should also render the case for the new start "
            "codon.",
            "Used NG_012337.1 instead of AB026906.1." "To be implemented.",
        ],
        # "input": "AB026906.1:c.1_4delinsTTGA", # yields TTG start codon
        "input": "NG_012337.1(NM_003002.2):c.1_4delinsTTGA",  # yields TTG start codon
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_003002.2):c.[1A>T;4G>A]",
                "NG_012337.1(NP_002993.1):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: start_codon_alt_start_yield_start",
            "Mutation creating new start codon should yield a p.? description, "
            "also with an alternative start codon. The visualisation should "
            "also render the new start codon.",
            "To be implemented.",
        ],
        "input": "NM_024426.4:c.1_4delinsATGA",  # yields ATG start codon
        "coding_protein_descriptions": {
            (
                "NM_024426.4(NM_024426.4):c.[1C>A;4C>A]",
                "NM_024426.4(NP_077744.3):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: protein_ext_stop",
            "Variant in stop codon where an alternative stop codon is found "
            "downstream in the RNA should yield `ext*P` where P is a position.",
        ],
        "input": "NM_000143.3:c.1531T>G",
        "coding_protein_descriptions": {
            (
                "NM_000143.3(NM_000143.3):c.1531T>G",
                "NM_000143.3(NP_000134.2):p.(*511Glyext*3)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.7125G>T",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_003002.2):c.274G>T",
                "NG_012337.1(NP_002993.1):p.(Asp92Tyr)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "LRG_24:g.5526_5533del",
        "coding_protein_descriptions": {
            (
                "LRG_24(t1):c.127_134del",
                "LRG_24(p1):p.(Gly43Argfs*65)",
            ),
            (
                "LRG_24(t2):c.127_134del",
                "LRG_24(p2):p.(Gly43Argfs*65)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_012459.2):c.4_5insGTA",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_012459.2):c.5_6insTAG",
                "NG_012337.1(NP_036591.2):p.(Arg2_Lys3insSer)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_012459.2):c.5_6delinsTAG",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_012459.2):c.5_6delinsTAG",
                "NG_012337.1(NP_036591.2):p.(Arg2Leufs*23)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_012459.2):c.4_6delinsGTA",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_012459.2):c.4_6delinsGTA",
                "NG_012337.1(NP_036591.2):p.(Arg2Val)",
            ),
        },
        "to_test": True,
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
