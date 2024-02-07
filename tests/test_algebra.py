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
                "supremal_lhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCCA",
                    "spdi": "LRG_24:5521:13:GCCCA",
                },
                "supremal_rhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCA",
                    "spdi": "LRG_24:5521:13:GCCA",
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
                "supremal_lhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCCA",
                    "spdi": "LRG_24:5521:13:GCCCA",
                },
                "supremal_rhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCA",
                    "spdi": "LRG_24:5521:13:GCCA",
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
                "supremal_lhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCCA",
                    "spdi": "LRG_24:5521:13:GCCCA",
                },
                "supremal_rhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCA",
                    "spdi": "LRG_24:5521:13:GCCA",
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
                "supremal_lhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCCA",
                    "spdi": "LRG_24:5521:13:GCCCA",
                },
                "supremal_rhs": {
                    "hgvs": "LRG_24:g.5522_5534delinsGCCA",
                    "spdi": "LRG_24:5521:13:GCCA",
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
                "supremal_lhs": {
                    "hgvs": "NG_008376.4:g.5519_5534delinsGAGTTATG",
                    "spdi": "NG_008376.4:5518:16:GAGTTATG",
                },
                "supremal_rhs": {
                    "hgvs": "NG_008376.4:g.5519_5534delinsGAGTTATG",
                    "spdi": "NG_008376.4:5518:16:GAGTTATG",
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
                "supremal_lhs": {
                    "hgvs": "NG_012337.3:g.273_275delinsAAA",
                    "spdi": "NG_012337.3:272:3:AAA",
                },
                "supremal_rhs": {
                    "hgvs": "NG_012337.3:g.274del",
                    "spdi": "NG_012337.3:273:1:",
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
                "supremal_lhs": {
                    "hgvs": "NG_012337.3:g.273_275delinsAAA",
                    "spdi": "NG_012337.3:272:3:AAA",
                },
                "supremal_rhs": {
                    "hgvs": "NG_012337.3:g.274del",
                    "spdi": "NG_012337.3:273:1:",
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
                "supremal_lhs": {
                    "hgvs": "NG_012337.3:g.273_275delinsAAA",
                    "spdi": "NG_012337.3:272:3:AAA",
                },
                "supremal_rhs": {
                    "hgvs": "NG_012337.3:g.274del",
                    "spdi": "NG_012337.3:273:1:",
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
                "supremal_lhs": {
                    "hgvs": "NG_012337.3:g.273_275delinsAAA",
                    "spdi": "NG_012337.3:272:3:AAA",
                },
                "supremal_rhs": {
                    "hgvs": "NG_012337.3:g.270_271del",
                    "spdi": "NG_012337.3:269:2:",
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
        (
            (
                "AAAAA",
                "sequence",
                "2_3insT",
                "variant",
                "ATAAAAA",
                "variant",
            ),
            {
                "errors": {
                    "rhs": [
                        {
                            "code": "ESYNTAXUC",
                            "details": "Unexpected character.",
                            "line": 1,
                            "column": 1,
                            "pos_in_stream": 0,
                            "unexpected_character": "A",
                            "description": "ATAAAAA",
                            "expecting": [
                                "'(=)' for predicted no changes",
                                "a number (to indicate a location or a length)",
                                "'*' or '-' for an outside CDS location",
                                "'[('",
                                "'(' for an uncertainty start or before a selector ID",
                                "?",
                                "'[' for multiple variants, insertions, or repeats",
                                "'pter' or 'qter'",
                                "'=' to indicate no changes",
                                "'(['",
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
    assert list(output["errors"].keys()) == list(expected["errors"].keys())
    assert (
        output["errors"]["rhs"][0]["pos_in_stream"]
        == expected["errors"]["rhs"][0]["pos_in_stream"]
    )
    assert (
        output["errors"]["rhs"][0]["unexpected_character"]
        == expected["errors"]["rhs"][0]["unexpected_character"]
    )
    assert (
        output["errors"]["rhs"][0]["description"]
        == expected["errors"]["rhs"][0]["description"]
    )
    assert set(output["errors"]["rhs"][0]["expecting"]) == set(
        expected["errors"]["rhs"][0]["expecting"]
    )
