import pytest

from mutalyzer.viewer import view_variants

from .commons import monkey_patches


@pytest.mark.parametrize(
    "description, output",
    [
        (
            "NM_003002.2:r.274g>u",
            {
                "views": [
                    {
                        "start": 0,
                        "end": 334,
                        "type": "outside",
                        "left": "gugggaauugucgcc",
                        "right": "ccuugcucugcgaug",
                    },
                    {
                        "description": "274g>u",
                        "start": 334,
                        "end": 335,
                        "type": "variant",
                        "deleted": {"sequence": "g"},
                        "inserted": {"sequence": "u", "length": 1},
                    },
                    {
                        "start": 335,
                        "end": 1382,
                        "type": "outside",
                        "left": "acuauucccuggcug",
                        "right": "aaaaaaaaaaaaaaa",
                    },
                ],
                "seq_length": 1382,
            },
        ),
        (
            "NG_012337.1(NM_012459.2):c.297_*1del",
            {
                "seq_length": 15948,
                "inverted": True,
                "views": [
                    {
                        "start": 0,
                        "end": 12499,
                        "type": "outside",
                        "left": "AGGCAGGTAAAATTT",
                        "right": "GAAAGGAGGGCAGTA",
                    },
                    {
                        "description": "297_*1del",
                        "start": 12499,
                        "end": 12501,
                        "type": "variant",
                        "deleted": {"sequence": "GG"},
                    },
                    {
                        "start": 12501,
                        "end": 15948,
                        "type": "outside",
                        "left": "CCATCCCCCAGGAGA",
                        "right": "GGTAGAACCAAGCCC",
                    },
                ],
            },
        ),
        (
            "NG_012337.1(NM_012459.2):c.[10del;23_25del;37_38insATA]",
            {
                "seq_length": 15948,
                "inverted": True,
                "views": [
                    {
                        "start": 0,
                        "end": 11035,
                        "type": "outside",
                        "left": "AGGCAGGTAAAATTT",
                        "right": "TCGCGCATGCGCAAA",
                    },
                    {
                        "description": "10del",
                        "start": 11035,
                        "end": 11036,
                        "type": "variant",
                        "deleted": {"sequence": "C"},
                    },
                    {
                        "start": 11036,
                        "end": 11048,
                        "type": "outside",
                        "sequence": "ACAGCTGTCGGA",
                    },
                    {
                        "description": "23_25del",
                        "start": 11048,
                        "end": 11051,
                        "type": "variant",
                        "deleted": {"sequence": "AGG"},
                    },
                    {
                        "start": 11051,
                        "end": 11063,
                        "type": "outside",
                        "sequence": "TGGCGAGCCTGA",
                    },
                    {
                        "description": "37_38insATA",
                        "start": 11063,
                        "end": 11063,
                        "type": "variant",
                        "inserted": {"sequence": "ATA", "length": 3},
                    },
                    {
                        "start": 11063,
                        "end": 15948,
                        "type": "outside",
                        "left": "GGCGAACAATGGCGG",
                        "right": "GGTAGAACCAAGCCC",
                    },
                ],
            },
        ),
        (
            "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGlu)",
            {
                "seq_length": 133,
                "views": [
                    {
                        "start": 0,
                        "end": 67,
                        "type": "outside",
                        "left": "MVRRFLVTLRIRRAC",
                        "right": "QRLGQQPLPRRPGHD",
                    },
                    {
                        "description": "Asp68_Gly69delinsGlu",
                        "start": 67,
                        "end": 69,
                        "type": "variant",
                        "deleted": {"sequence": "DG"},
                        "inserted": {"sequence": "E", "length": 1},
                    },
                    {
                        "start": 69,
                        "end": 133,
                        "type": "outside",
                        "left": "QRPSGGAAAAPRRGA",
                        "right": "GRARCLGPSARGPG*",
                    },
                ],
            },
        ),
        (
            "NG_007485.1(NP_000068.1):p.(Met54delinsIleSer)",
            {
                "seq_length": 157,
                "views": [
                    {
                        "start": 0,
                        "end": 53,
                        "type": "outside",
                        "left": "MEPAAGSSMEPSADW",
                        "right": "NAPNSYGRRPIQVMM",
                    },
                    {
                        "description": "Met54delinsIleSer",
                        "start": 53,
                        "end": 54,
                        "type": "variant",
                        "deleted": {"sequence": "M"},
                        "inserted": {"sequence": "IS", "length": 2},
                    },
                    {
                        "start": 54,
                        "end": 157,
                        "type": "outside",
                        "left": "GSARVAELLLLHGAE",
                        "right": "ARIDAAEGPSDIPD*",
                    },
                ],
            },
        ),
    ],
)
def test_view_variants(description, output):
    assert view_variants(description) == output


@pytest.mark.parametrize(
    "description, errors",
    [
        (
            "NG_007485.1(NP_000068.1):p.(Met54Ilefs*66)",
            [
                {
                    "code": "EVARIANTNOTSUPPORTED",
                    "details": "Variant Met54[Ile;66] type frame_shift not supported.",
                    "paths": [["variants", 0]],
                }
            ],
        )
    ],
)
def test_view_variants_errors(description, errors):
    assert view_variants(description)["errors"] == errors
