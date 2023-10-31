import pytest

from mutalyzer.algebra import _get_hgvs_and_variant, _get_id, compare
from mutalyzer.reference import retrieve_reference

from .commons import code_in, monkey_patches


@pytest.mark.parametrize(
    "reference_id",
    ["NM_012459.2", "LRG_24"],
)
def test_get_id(reference_id):
    assert _get_id(reference_id) == {
        "input": reference_id,
        "type": "id",
        "reference": {"id": reference_id},
        "annotations": {"id": reference_id},
        "sequence": retrieve_reference(reference_id)[0]["sequence"]["seq"],
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
                "view": {
                    "views": [
                        {
                            "start": 0,
                            "end": 30,
                            "type": "outside",
                            "sequence": "AAGTCGAGAGGCGGTGCACACCCGTCGCGC",
                        },
                        {
                            "description": "1del",
                            "start": 30,
                            "end": 31,
                            "type": "variant",
                            "deleted": {"sequence": "A"},
                        },
                        {
                            "start": 31,
                            "end": 823,
                            "type": "outside",
                            "left": "TGCGCAAACACAGCT",
                            "right": "AAAAAAAAAAAAAAA",
                        },
                    ],
                    "seq_length": 823,
                },
            },
        )
    ],
)
def test_get_hgvs(description, expected):
    assert _get_hgvs_and_variant(description) == expected


@pytest.mark.parametrize(
    "params, expected",
    [
        (  # 0
            ("AAAAA", "sequence", "ATAAAAA", "sequence", "AATAAA", "sequence"),
            {
                "relation": "disjoint",
                "influence_lhs": {"min_pos": 0, "max_pos": 5},
                "influence_rhs": {"min_pos": 2, "max_pos": 2},
                "supremal_lhs": {"hgvs": "1_5delinsATAAAAA", "spdi": ":0:5:ATAAAAA"},
                "supremal_rhs": {"hgvs": "2_3insT", "spdi": ":2:0:T"},
                "view_lhs_supremal": {
                    "seq_length": 5,
                    "views": [
                        {"start": 0, "end": 0, "type": "outside"},
                        {
                            "description": "1_5delinsATAAAAA",
                            "start": 0,
                            "end": 5,
                            "type": "variant",
                            "deleted": {"sequence": "AAAAA"},
                            "inserted": {"sequence": "ATAAAAA", "length": 7},
                        },
                        {"start": 5, "end": 5, "type": "outside"},
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 5,
                    "views": [
                        {"start": 0, "end": 2, "type": "outside", "sequence": "AA"},
                        {
                            "description": "2_3insT",
                            "start": 2,
                            "end": 2,
                            "type": "variant",
                            "inserted": {"sequence": "T", "length": 1},
                        },
                        {"start": 2, "end": 5, "type": "outside", "sequence": "AAA"},
                    ],
                },
            },
        ),
        (  # 1
            ("AAAAA", "sequence", "ATAAAAA", "sequence", "2_3insT", "variant"),
            {
                "relation": "disjoint",
                "influence_lhs": {"min_pos": 0, "max_pos": 5},
                "influence_rhs": {"min_pos": 2, "max_pos": 2},
                "supremal_lhs": {"hgvs": "1_5delinsATAAAAA", "spdi": ":0:5:ATAAAAA"},
                "supremal_rhs": {"hgvs": "2_3insT", "spdi": ":2:0:T"},
                "view_rhs": {
                    "views": [
                        {"start": 0, "end": 2, "type": "outside", "sequence": "AA"},
                        {
                            "description": "2_3insT",
                            "start": 2,
                            "end": 2,
                            "type": "variant",
                            "inserted": {"sequence": "T", "length": 1},
                        },
                        {"start": 2, "end": 5, "type": "outside", "sequence": "AAA"},
                    ],
                    "seq_length": 5,
                },
                "view_lhs_supremal": {
                    "seq_length": 5,
                    "views": [
                        {"start": 0, "end": 0, "type": "outside"},
                        {
                            "description": "1_5delinsATAAAAA",
                            "start": 0,
                            "end": 5,
                            "type": "variant",
                            "deleted": {"sequence": "AAAAA"},
                            "inserted": {"sequence": "ATAAAAA", "length": 7},
                        },
                        {"start": 5, "end": 5, "type": "outside"},
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 5,
                    "views": [
                        {"start": 0, "end": 2, "type": "outside", "sequence": "AA"},
                        {
                            "description": "2_3insT",
                            "start": 2,
                            "end": 2,
                            "type": "variant",
                            "inserted": {"sequence": "T", "length": 1},
                        },
                        {"start": 2, "end": 5, "type": "outside", "sequence": "AAA"},
                    ],
                },
            },
        ),
        (  # 2
            (
                None,
                None,
                "LRG_24:g.5525_5532del",
                "hgvs",
                "LRG_24:g.5525_5533del",
                "hgvs",
            ),
            {
                "relation": "is_contained",
                "influence_lhs": {"min_pos": 5521, "max_pos": 5534},
                "influence_rhs": {"min_pos": 5521, "max_pos": 5534},
                "supremal_lhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCCA",
                    "spdi": "LRG_24:5521:13:GCCCA",
                },
                "supremal_rhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCA",
                    "spdi": "LRG_24:5521:13:GCCA",
                },
                "view_lhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "GATGCCCGGCCTGCC",
                        },
                        {
                            "description": "5525_5532del",
                            "start": 5524,
                            "end": 5532,
                            "type": "variant",
                            "deleted": {"sequence": "CGGGGCAC"},
                        },
                        {
                            "start": 5532,
                            "end": 11486,
                            "type": "outside",
                            "left": "CAGGGAAGGATGGGT",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                    "seq_length": 11486,
                },
                "view_rhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "GATGCCCGGCCTGCC",
                        },
                        {
                            "description": "5525_5533del",
                            "start": 5524,
                            "end": 5533,
                            "type": "variant",
                            "deleted": {"sequence": "CGGGGCACC"},
                        },
                        {
                            "start": 5533,
                            "end": 11486,
                            "type": "outside",
                            "left": "AGGGAAGGATGGGTA",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                    "seq_length": 11486,
                },
                "view_lhs_supremal": {
                    "seq_length": 11486,
                    "views": [
                        {
                            "start": 0,
                            "end": 5521,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "AGGGATGCCCGGCCT",
                        },
                        {
                            "description": "5522_5534delinsGCCCA",
                            "start": 5521,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GCCCGGGGCACCA"},
                            "inserted": {"sequence": "GCCCA", "length": 5},
                        },
                        {
                            "start": 5534,
                            "end": 11486,
                            "type": "outside",
                            "left": "GGGAAGGATGGGTAC",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 11486,
                    "views": [
                        {
                            "start": 0,
                            "end": 5521,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "AGGGATGCCCGGCCT",
                        },
                        {
                            "description": "5522_5534delinsGCCA",
                            "start": 5521,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GCCCGGGGCACCA"},
                            "inserted": {"sequence": "GCCA", "length": 4},
                        },
                        {
                            "start": 5534,
                            "end": 11486,
                            "type": "outside",
                            "left": "GGGAAGGATGGGTAC",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                },
            },
        ),
        (  # 3
            (
                None,
                None,
                "LRG_24:g.5525_5532del",
                "hgvs",
                "LRG_24:g.5525_5533del",
                "hgvs",
            ),
            {
                "relation": "is_contained",
                "influence_lhs": {"min_pos": 5521, "max_pos": 5534},
                "influence_rhs": {"min_pos": 5521, "max_pos": 5534},
                "supremal_lhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCCA",
                    "spdi": "LRG_24:5521:13:GCCCA",
                },
                "supremal_rhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCA",
                    "spdi": "LRG_24:5521:13:GCCA",
                },
                "view_lhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "GATGCCCGGCCTGCC",
                        },
                        {
                            "description": "5525_5532del",
                            "start": 5524,
                            "end": 5532,
                            "type": "variant",
                            "deleted": {"sequence": "CGGGGCAC"},
                        },
                        {
                            "start": 5532,
                            "end": 11486,
                            "type": "outside",
                            "left": "CAGGGAAGGATGGGT",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                    "seq_length": 11486,
                },
                "view_rhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "GATGCCCGGCCTGCC",
                        },
                        {
                            "description": "5525_5533del",
                            "start": 5524,
                            "end": 5533,
                            "type": "variant",
                            "deleted": {"sequence": "CGGGGCACC"},
                        },
                        {
                            "start": 5533,
                            "end": 11486,
                            "type": "outside",
                            "left": "AGGGAAGGATGGGTA",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                    "seq_length": 11486,
                },
                "view_lhs_supremal": {
                    "seq_length": 11486,
                    "views": [
                        {
                            "start": 0,
                            "end": 5521,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "AGGGATGCCCGGCCT",
                        },
                        {
                            "description": "5522_5534delinsGCCCA",
                            "start": 5521,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GCCCGGGGCACCA"},
                            "inserted": {"sequence": "GCCCA", "length": 5},
                        },
                        {
                            "start": 5534,
                            "end": 11486,
                            "type": "outside",
                            "left": "GGGAAGGATGGGTAC",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 11486,
                    "views": [
                        {
                            "start": 0,
                            "end": 5521,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "AGGGATGCCCGGCCT",
                        },
                        {
                            "description": "5522_5534delinsGCCA",
                            "start": 5521,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GCCCGGGGCACCA"},
                            "inserted": {"sequence": "GCCA", "length": 4},
                        },
                        {
                            "start": 5534,
                            "end": 11486,
                            "type": "outside",
                            "left": "GGGAAGGATGGGTAC",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                },
            },
        ),
        (  # 4
            (
                "LRG_24",
                "id",
                "LRG_24:g.5525_5532del",
                "hgvs",
                "LRG_24:g.5525_5533del",
                "hgvs",
            ),
            {
                "relation": "is_contained",
                "influence_lhs": {"min_pos": 5521, "max_pos": 5534},
                "influence_rhs": {"min_pos": 5521, "max_pos": 5534},
                "supremal_lhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCCA",
                    "spdi": "LRG_24:5521:13:GCCCA",
                },
                "supremal_rhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCA",
                    "spdi": "LRG_24:5521:13:GCCA",
                },
                "view_lhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "GATGCCCGGCCTGCC",
                        },
                        {
                            "description": "5525_5532del",
                            "start": 5524,
                            "end": 5532,
                            "type": "variant",
                            "deleted": {"sequence": "CGGGGCAC"},
                        },
                        {
                            "start": 5532,
                            "end": 11486,
                            "type": "outside",
                            "left": "CAGGGAAGGATGGGT",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                    "seq_length": 11486,
                },
                "view_rhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "GATGCCCGGCCTGCC",
                        },
                        {
                            "description": "5525_5533del",
                            "start": 5524,
                            "end": 5533,
                            "type": "variant",
                            "deleted": {"sequence": "CGGGGCACC"},
                        },
                        {
                            "start": 5533,
                            "end": 11486,
                            "type": "outside",
                            "left": "AGGGAAGGATGGGTA",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                    "seq_length": 11486,
                },
                "view_lhs_supremal": {
                    "seq_length": 11486,
                    "views": [
                        {
                            "start": 0,
                            "end": 5521,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "AGGGATGCCCGGCCT",
                        },
                        {
                            "description": "5522_5534delinsGCCCA",
                            "start": 5521,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GCCCGGGGCACCA"},
                            "inserted": {"sequence": "GCCCA", "length": 5},
                        },
                        {
                            "start": 5534,
                            "end": 11486,
                            "type": "outside",
                            "left": "GGGAAGGATGGGTAC",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 11486,
                    "views": [
                        {
                            "start": 0,
                            "end": 5521,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "AGGGATGCCCGGCCT",
                        },
                        {
                            "description": "5522_5534delinsGCCA",
                            "start": 5521,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GCCCGGGGCACCA"},
                            "inserted": {"sequence": "GCCA", "length": 4},
                        },
                        {
                            "start": 5534,
                            "end": 11486,
                            "type": "outside",
                            "left": "GGGAAGGATGGGTAC",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                },
            },
        ),
        (  # 5
            (
                "LRG_24",
                "id",
                "LRG_24:g.5525_5532del",
                "hgvs",
                "LRG_24:g.5525_5533del",
                "hgvs",
            ),
            {
                "relation": "is_contained",
                "influence_lhs": {"min_pos": 5521, "max_pos": 5534},
                "influence_rhs": {"min_pos": 5521, "max_pos": 5534},
                "supremal_lhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCCA",
                    "spdi": "LRG_24:5521:13:GCCCA",
                },
                "supremal_rhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCA",
                    "spdi": "LRG_24:5521:13:GCCA",
                },
                "view_lhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "GATGCCCGGCCTGCC",
                        },
                        {
                            "description": "5525_5532del",
                            "start": 5524,
                            "end": 5532,
                            "type": "variant",
                            "deleted": {"sequence": "CGGGGCAC"},
                        },
                        {
                            "start": 5532,
                            "end": 11486,
                            "type": "outside",
                            "left": "CAGGGAAGGATGGGT",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                    "seq_length": 11486,
                },
                "view_rhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "GATGCCCGGCCTGCC",
                        },
                        {
                            "description": "5525_5533del",
                            "start": 5524,
                            "end": 5533,
                            "type": "variant",
                            "deleted": {"sequence": "CGGGGCACC"},
                        },
                        {
                            "start": 5533,
                            "end": 11486,
                            "type": "outside",
                            "left": "AGGGAAGGATGGGTA",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                    "seq_length": 11486,
                },
                "view_lhs_supremal": {
                    "seq_length": 11486,
                    "views": [
                        {
                            "start": 0,
                            "end": 5521,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "AGGGATGCCCGGCCT",
                        },
                        {
                            "description": "5522_5534delinsGCCCA",
                            "start": 5521,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GCCCGGGGCACCA"},
                            "inserted": {"sequence": "GCCCA", "length": 5},
                        },
                        {
                            "start": 5534,
                            "end": 11486,
                            "type": "outside",
                            "left": "GGGAAGGATGGGTAC",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 11486,
                    "views": [
                        {
                            "start": 0,
                            "end": 5521,
                            "type": "outside",
                            "left": "GTTCACACTTCTCTC",
                            "right": "AGGGATGCCCGGCCT",
                        },
                        {
                            "description": "5522_5534delinsGCCA",
                            "start": 5521,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GCCCGGGGCACCA"},
                            "inserted": {"sequence": "GCCA", "length": 4},
                        },
                        {
                            "start": 5534,
                            "end": 11486,
                            "type": "outside",
                            "left": "GGGAAGGATGGGTAC",
                            "right": "CATGTATACACATAC",
                        },
                    ],
                },
            },
        ),
        (  # 6
            (
                None,
                None,
                "NG_008376.4:g.5525_5532del",
                "hgvs",
                "LRG_303:g.5525_5532del",
                "hgvs",
            ),
            {
                "relation": "equivalent",
                "influence_lhs": {"min_pos": 5518, "max_pos": 5534},
                "influence_rhs": {"min_pos": 5518, "max_pos": 5534},
                "supremal_lhs": {
                    "hgvs": "NG_008376.4:g.5519_5534delinsGAGTTATG",
                    "spdi": "NG_008376.4:5518:16:GAGTTATG",
                },
                "supremal_rhs": {
                    "hgvs": "NG_008376.4:g.5519_5534delinsGAGTTATG",
                    "spdi": "NG_008376.4:5518:16:GAGTTATG",
                },
                "view_lhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTGTCTGCCAAGGGT",
                            "right": "ATGGAAGATGAGTTA",
                        },
                        {
                            "description": "5525_5532del",
                            "start": 5524,
                            "end": 5532,
                            "type": "variant",
                            "deleted": {"sequence": "GTCCTGAG"},
                        },
                        {
                            "start": 5532,
                            "end": 11312,
                            "type": "outside",
                            "left": "TGCCGTTTAAATCAC",
                            "right": "GGAAGAAGAGGGATC",
                        },
                    ],
                    "seq_length": 11312,
                },
                "view_rhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 5524,
                            "type": "outside",
                            "left": "GTGTCTGCCAAGGGT",
                            "right": "ATGGAAGATGAGTTA",
                        },
                        {
                            "description": "5525_5532del",
                            "start": 5524,
                            "end": 5532,
                            "type": "variant",
                            "deleted": {"sequence": "GTCCTGAG"},
                        },
                        {
                            "start": 5532,
                            "end": 11312,
                            "type": "outside",
                            "left": "TGCCGTTTAAATCAC",
                            "right": "GGAAGAAGAGGGATC",
                        },
                    ],
                    "seq_length": 11312,
                },
                "view_lhs_supremal": {
                    "seq_length": 11312,
                    "views": [
                        {
                            "start": 0,
                            "end": 5518,
                            "type": "outside",
                            "left": "GTGTCTGCCAAGGGT",
                            "right": "TTGGGAATGGAAGAT",
                        },
                        {
                            "description": "5519_5534delinsGAGTTATG",
                            "start": 5518,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GAGTTAGTCCTGAGTG"},
                            "inserted": {"sequence": "GAGTTATG", "length": 8},
                        },
                        {
                            "start": 5534,
                            "end": 11312,
                            "type": "outside",
                            "left": "CCGTTTAAATCACGA",
                            "right": "GGAAGAAGAGGGATC",
                        },
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 11312,
                    "views": [
                        {
                            "start": 0,
                            "end": 5518,
                            "type": "outside",
                            "left": "GTGTCTGCCAAGGGT",
                            "right": "TTGGGAATGGAAGAT",
                        },
                        {
                            "description": "5519_5534delinsGAGTTATG",
                            "start": 5518,
                            "end": 5534,
                            "type": "variant",
                            "deleted": {"sequence": "GAGTTAGTCCTGAGTG"},
                            "inserted": {"sequence": "GAGTTATG", "length": 8},
                        },
                        {
                            "start": 5534,
                            "end": 11312,
                            "type": "outside",
                            "left": "CCGTTTAAATCACGA",
                            "right": "GGAAGAAGAGGGATC",
                        },
                    ],
                },
            },
        ),
        (  # 7
            (
                None,
                None,
                "NG_012337.3:g.274>A",
                "hgvs",
                "274del",
                "variant",
            ),
            {
                "relation": "contains",
                "influence_lhs": {"min_pos": 272, "max_pos": 275},
                "influence_rhs": {"min_pos": 273, "max_pos": 274},
                "supremal_lhs": {
                    "hgvs": "NG_012337.3:g.273_275delinsAAA",
                    "spdi": "NG_012337.3:272:3:AAA",
                },
                "supremal_rhs": {
                    "hgvs": "NG_012337.3:g.274del",
                    "spdi": "NG_012337.3:273:1:",
                },
                "view_lhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274>A",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                            "inserted": {"sequence": "A", "length": 1},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                    "seq_length": 39784,
                },
                "view_rhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274del",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                    "seq_length": 39784,
                },
                "view_lhs_supremal": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 272,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "ACTTTTACGAAGAAT",
                        },
                        {
                            "description": "273_275delinsAAA",
                            "start": 272,
                            "end": 275,
                            "type": "variant",
                            "deleted": {"sequence": "ATA"},
                            "inserted": {"sequence": "AAA", "length": 3},
                        },
                        {
                            "start": 275,
                            "end": 39784,
                            "type": "outside",
                            "left": "CTTGCCATCAAAAAT",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274del",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
            },
        ),
        (  # 8
            (
                None,
                None,
                "NG_012337.3:g.274>A",
                "hgvs",
                "274del",
                "variant",
            ),
            {
                "relation": "contains",
                "influence_lhs": {"min_pos": 272, "max_pos": 275},
                "influence_rhs": {"min_pos": 273, "max_pos": 274},
                "supremal_lhs": {
                    "hgvs": "NG_012337.3:g.273_275delinsAAA",
                    "spdi": "NG_012337.3:272:3:AAA",
                },
                "supremal_rhs": {
                    "hgvs": "NG_012337.3:g.274del",
                    "spdi": "NG_012337.3:273:1:",
                },
                "view_lhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274>A",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                            "inserted": {"sequence": "A", "length": 1},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                    "seq_length": 39784,
                },
                "view_rhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274del",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                    "seq_length": 39784,
                },
                "view_lhs_supremal": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 272,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "ACTTTTACGAAGAAT",
                        },
                        {
                            "description": "273_275delinsAAA",
                            "start": 272,
                            "end": 275,
                            "type": "variant",
                            "deleted": {"sequence": "ATA"},
                            "inserted": {"sequence": "AAA", "length": 3},
                        },
                        {
                            "start": 275,
                            "end": 39784,
                            "type": "outside",
                            "left": "CTTGCCATCAAAAAT",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274del",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
            },
        ),
        (  # 9
            (
                "NG_012337.3",
                "id",
                "274>A",
                "variant",
                "274del",
                "variant",
            ),
            {
                "relation": "contains",
                "influence_lhs": {"min_pos": 272, "max_pos": 275},
                "influence_rhs": {"min_pos": 273, "max_pos": 274},
                "supremal_lhs": {
                    "hgvs": "NG_012337.3:g.273_275delinsAAA",
                    "spdi": "NG_012337.3:272:3:AAA",
                },
                "supremal_rhs": {
                    "hgvs": "NG_012337.3:g.274del",
                    "spdi": "NG_012337.3:273:1:",
                },
                "view_lhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274>A",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                            "inserted": {"sequence": "A", "length": 1},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                    "seq_length": 39784,
                },
                "view_rhs": {
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274del",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                    "seq_length": 39784,
                },
                "view_lhs_supremal": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 272,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "ACTTTTACGAAGAAT",
                        },
                        {
                            "description": "273_275delinsAAA",
                            "start": 272,
                            "end": 275,
                            "type": "variant",
                            "deleted": {"sequence": "ATA"},
                            "inserted": {"sequence": "AAA", "length": 3},
                        },
                        {
                            "start": 275,
                            "end": 39784,
                            "type": "outside",
                            "left": "CTTGCCATCAAAAAT",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274del",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
            },
        ),
        (  # 10
            (
                "NG_012337.3",
                "id",
                "274>A",
                "variant",
                "270_271del",
                "variant",
            ),
            {
                "relation": "disjoint",
                "influence_lhs": {"min_pos": 272, "max_pos": 275},
                "influence_rhs": {"min_pos": 269, "max_pos": 271},
                "supremal_lhs": {
                    "hgvs": "NG_012337.3:g.273_275delinsAAA",
                    "spdi": "NG_012337.3:272:3:AAA",
                },
                "supremal_rhs": {
                    "hgvs": "NG_012337.3:g.270_271del",
                    "spdi": "NG_012337.3:269:2:",
                },
                "view_lhs": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 273,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "CTTTTACGAAGAATA",
                        },
                        {
                            "description": "274>A",
                            "start": 273,
                            "end": 274,
                            "type": "variant",
                            "deleted": {"sequence": "T"},
                            "inserted": {"sequence": "A", "length": 1},
                        },
                        {
                            "start": 274,
                            "end": 39784,
                            "type": "outside",
                            "left": "ACTTGCCATCAAAAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
                "view_rhs": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 269,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "TTAACTTTTACGAAG",
                        },
                        {
                            "description": "270_271del",
                            "start": 269,
                            "end": 271,
                            "type": "variant",
                            "deleted": {"sequence": "AA"},
                        },
                        {
                            "start": 271,
                            "end": 39784,
                            "type": "outside",
                            "left": "TATACTTGCCATCAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
                "view_lhs_supremal": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 272,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "ACTTTTACGAAGAAT",
                        },
                        {
                            "description": "273_275delinsAAA",
                            "start": 272,
                            "end": 275,
                            "type": "variant",
                            "deleted": {"sequence": "ATA"},
                            "inserted": {"sequence": "AAA", "length": 3},
                        },
                        {
                            "start": 275,
                            "end": 39784,
                            "type": "outside",
                            "left": "CTTGCCATCAAAAAT",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
                "view_rhs_supremal": {
                    "seq_length": 39784,
                    "views": [
                        {
                            "start": 0,
                            "end": 269,
                            "type": "outside",
                            "left": "GGGCTTGGTTCTACC",
                            "right": "TTAACTTTTACGAAG",
                        },
                        {
                            "description": "270_271del",
                            "start": 269,
                            "end": 271,
                            "type": "variant",
                            "deleted": {"sequence": "AA"},
                        },
                        {
                            "start": 271,
                            "end": 39784,
                            "type": "outside",
                            "left": "TATACTTGCCATCAA",
                            "right": "AAATTACTCAAGGAA",
                        },
                    ],
                },
            },
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
                            "details": "hgvs not valid.",
                            "options": ["sequence", "id"],
                        }
                    ],
                    "lhs_type": [
                        {
                            "code": "EINVALIDINPUT",
                            "details": "HGVS not valid.",
                            "options": ["sequence", "variant", "hgvs"],
                        }
                    ],
                    "rhs_type": [
                        {
                            "code": "EINVALIDINPUT",
                            "details": "varianT not valid.",
                            "options": ["sequence", "variant", "hgvs"],
                        }
                    ],
                }
            },
        ),
    ],
)
def test_compare_errors(params, expected):
    assert compare(*params) == expected


@pytest.mark.parametrize(
    "params, expected",
    [
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
                            "pos_in_stream": 6,
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
def test_compare_errors_syntax(params, expected):
    output = compare(*params)
    assert list(output.keys()) == ["errors"]
    assert list(output["errors"].keys()) == ["rhs"]
    assert output["errors"]["rhs"][0]["pos_in_stream"] == 6
    assert output["errors"]["rhs"][0]["unexpected_character"] == "e"
    assert output["errors"]["rhs"][0]["description"] == "1delete"
    assert set(output["errors"]["rhs"][0]["expecting"]) == {
        "'(' for an uncertainty start or before a selector ID",
        "':' between the reference part and the coordinate system",
        "a reference / selector ID",
    }
