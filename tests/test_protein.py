import pytest

from mutalyzer.normalizer import normalize

from .commons import code_in, monkey_patches

TESTS = [
    {
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.=",
        "normalized": "NG_012337.1(NP_002993.1):p.=",
        "to_test": True,
    },
    {
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.(=)",
        "normalized": "NG_012337.1(NP_002993.1):p.(=)",
        "to_test": True,
    },
    {
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.=",
        "normalized": "NG_012337.1(NP_002993.1):p.=",
        "to_test": True,
    },
    {
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92=",
        "normalized": "NG_012337.1(NP_002993.1):p.Asp92=",
        "to_test": True,
    },
    {
        "keywords": ["protein", "equal", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92Asp",
        "normalized": "NG_012337.1(NP_002993.1):p.=",
        "to_test": True,
    },
    {
        "keywords": ["protein", "substitution", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92Tyr",
        "normalized": "NG_012337.1(NP_002993.1):p.Asp92Tyr",
        "to_test": True,
    },
    {
        "keywords": ["protein", "delins", "genomic"],
        "input": "NG_012337.1(NP_002993.1):p.Asp92delinsTyr",
        "normalized": "NG_012337.1(NP_002993.1):p.Asp92Tyr",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
            "delins",
            "substitution",
            "mRNA",
            "translation exception",
        ],
        "input": "NM_005410.4:p.U59delinsTyr",
        "normalized": "NM_005410.4(NP_005401.3):p.Sec59Tyr",
        "to_test": True,
    },
    {
        "keywords": ["protein", "delins", "mRNA", "translation exception", "error"],
        "input": "NM_005410.4:p.D59delinsA",
        "errors": [
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "Amino acid sequence `D` was not found in the reference sequence; `U` was found instead.",
                "paths": [["variants", 0, "location"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein", "delins", "protein", "translation exception", "error"],
        "input": "NP_005401.3:p.D59delinsA",
        "errors": [
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "Amino acid sequence `D` was not found in the reference sequence; `U` was found instead.",
                "paths": [["variants", 0, "location"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
            "delins",
            "mRNA",
            "reference insertion",
            "translation exception",
        ],
        "input": "NM_005410.4:p.U59delinsNP_002993.1:H50_53",
        "normalized": "NM_005410.4(NP_005401.3):p.Sec59delinsHisLeuSerPro",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
            "delins",
            "mRNA",
            "reference insertion",
            "translation exception",
            "error",
        ],
        "input": "NM_005410.4:p.D59delinsNP_002993.1:Asp50_53",
        "errors": [
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "Amino acid sequence `D` was not found in the reference sequence; `U` was found instead.",
                "paths": [["variants", 0, "location"]],
            },
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "Amino acid sequence `D` was not found in the reference sequence; `H` was found instead.",
                "paths": [["variants", 0, "inserted", 0, "location", "start"]],
            },
        ],
        "to_test": True,
    },
    {
        "keywords": [
            "protein missing mrna for backtranslation",
        ],
        "input": "YP_009725300.1:p.(Leu360Ter)",
        "normalized": "YP_009725300.1:p.(Leu360Ter)",
        "to_test": True,
    },
    {
        "keywords": [
            "protein reverse strand",
        ],
        "input": "NG_008835.1(NP_071436.1):p.(Lys34Val)",
        "normalized": "NG_008835.1(NP_071436.1):p.(Lys34Val)",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
        ],
        "input": "NP_071436.1:p.(Lys34Val)",
        "normalized": "NP_071436.1:p.(Lys34Val)",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
        ],
        "input": "NP_071436.1:p.(Lys34Xaa)",
        "normalized": "NP_071436.1:p.(Lys34Xaa)",
        "to_test": True,
    },
    {
        "keywords": [
            "protein",
        ],
        "input": "NP_071436.1:p.(K34X)",
        "normalized": "NP_071436.1:p.(Lys34Xaa)",
        "to_test": True,
    },
    {
        "keywords": ["protein", "equal", "mitochondrial", "codon", "translation table"],
        # ENST00000361899.2 (MT-ATP6) contains an internal in frame TGA
        # codon, which is Trp under the vertebrate mitochondrial translation
        # table 2, but a premature stop under the standard table (1).
        "input": "ENST00000361899.2:p.Trp48=",
        "normalized": "ENST00000361899.2(ENSP00000354632):p.Trp48=",
        "to_test": True,
    },
    # NP_002993.1 has 159 amino acids, its stop codon is Ter160.
    {
        "keywords": ["protein", "stop codon"],
        "input": "NP_002993.1:p.Ter160=",
        "normalized": "NP_002993.1:p.Ter160=",
        "to_test": True,
    },
    {
        "keywords": ["protein", "stop codon", "extension"],
        "input": "NP_002993.1:p.Ter160Arg",
        "normalized": "NP_002993.1:p.Ter160ArgextTer?",
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.*160Glnext*17",
        "normalized": "NP_002993.1:p.Ter160GlnextTer17",
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.(Ter160GlnextTer?)",
        "normalized": "NP_002993.1:p.(Ter160GlnextTer?)",
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.Ter160GlnextTer1000",
        "normalized": "NP_002993.1:p.Ter160GlnextTer1000",
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.Met1ext-5",
        "normalized": "NP_002993.1:p.Met1ext-5",
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.Leu159GlnextTer17",
        "errors": [
            {
                "code": "EEXTENSION",
                "details": "Extension `Leu159GlnextTer17` is invalid; a C-terminal extension starts at the "
                "stop codon, not at `Leu159`.",
                "paths": [["variants", 0]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.Met1ext-0",
        "errors": [
            {
                "code": "EEXTENSION",
                "details": "Extension `Met1ext-0` is invalid; the new initiation position must be negative, "
                "not zero.",
                "paths": [["variants", 0]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.Ter160GlnGlyextTer5",
        "errors": [
            {
                "code": "EEXTENSION",
                "details": "Extension `Ter160GlnGlyextTer5` is invalid; the stop codon must become one "
                "amino acid, not `GlnGly`.",
                "paths": [["variants", 0]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.Ala1ext-5",
        "errors": [
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "Amino acid sequence `A` was not found in the reference sequence; `M` was found "
                "instead.",
                "paths": [["variants", 0, "location"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.Ter159GlnextTer17",
        "errors": [
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "Amino acid sequence `*` was not found in the reference sequence; `L` was found "
                "instead.",
                "paths": [["variants", 0, "location"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein", "extension"],
        "input": "NP_002993.1:p.Ter161GlnextTer17",
        "errors": [
            {
                "code": "EOUTOFBOUNDARY",
                "details": "Position `Ter161` is invalid; it is 1 amino acid after the end of the sequence "
                "(sequence length is 159 amino acids and the stop codon).",
                "paths": [["variants", 0, "location"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein"],
        "input": "NP_002993.1:p.Ala1Gly",
        "errors": [
            {
                "code": "EAMINOACIDMISMATCH",
                "details": "Amino acid sequence `A` was not found in the reference sequence; `M` was found "
                "instead.",
                "paths": [["variants", 0, "location"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein", "stop codon", "extension"],
        "input": "NP_002993.1:p.Leu159_Ter160insGly",
        "normalized": "NP_002993.1:p.Ter160GlyextTer1",
        "to_test": True,
    },
    {
        "keywords": ["protein", "stop codon"],
        "input": "NP_002993.1:p.Ter160del",
        "errors": [
            {
                "code": "ESTOPLOSS",
                "details": "Variant `Ter160del` removes the stop codon; the resulting extension "
                "cannot be described from the protein sequence.",
                "paths": [["variants"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein", "stop codon"],
        "input": "NP_002993.1:p.Ter161Gly",
        "errors": [
            {
                "code": "EOUTOFBOUNDARY",
                "details": "Position `Ter161` is invalid; it is 1 amino acid after the end of the sequence "
                "(sequence length is 159 amino acids and the stop codon).",
                "paths": [["variants", 0, "location"]],
            }
        ],
        "to_test": True,
    },
]


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test"):
            if test.get(t_type):
                output.append((test["input"], test[t_type]))
            else:
                output.append((test["input"], None))
    return output


@pytest.mark.parametrize(
    "input_description, normalized",
    get_tests(TESTS, "normalized"),
)
def test_protein(input_description, normalized):
    d = normalize(input_description)
    assert d.get("normalized_description") == normalized


@pytest.mark.parametrize(
    "input_description, errors",
    get_tests(TESTS, "errors"),
)
def test_protein_errors(input_description, errors):
    d = normalize(input_description)
    assert d.get("errors") == errors


@pytest.mark.parametrize(
    "variant, reason",
    [
        (
            "Leu159ext-5",
            "an N-terminal extension starts at `Met1`, not at `Leu159`",
        ),
        (
            "Met1_Leu159ext-5",
            "an N-terminal extension starts at `Met1`, not at `Met1_Leu159`",
        ),
        (
            "Ter160TerextTer5",
            "the stop codon must become one amino acid, not `Ter`",
        ),
        (
            "Ter160GlnextGly5",
            "the extension must end with `Ter`, not `Gly`",
        ),
        (
            "Ter160GlnextTer",
            "the extension length must be a positive number or `?`",
        ),
        (
            "Ter160GlnextTer0",
            "the extension length must be a positive number or `?`",
        ),
        (
            "Ter160GlnextTer17_20",
            "the extension length must be a positive number or `?`",
        ),
    ],
)
def test_extension_error_reason(variant, reason):
    result = normalize(f"NP_002993.1:p.{variant}")
    assert result.get("normalized_description") is None
    assert result["errors"] == [
        {
            "code": "EEXTENSION",
            "details": f"Extension `{variant}` is invalid; {reason}.",
            "paths": [["variants", 0]],
        }
    ]
