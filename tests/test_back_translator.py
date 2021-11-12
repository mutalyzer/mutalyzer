import pytest

from mutalyzer.back_translator import back_translate

from .commons import code_in, patch_retriever


@pytest.mark.parametrize(
    "protein_description, back_translated_descriptions",
    [
        (
            "NM_003002.2:p.Asp92Tyr",
            ["NM_003002.2:c.(274G>T)"],
        ),
        (
            "NG_012337.1(NP_002993.1):p.Asp92Tyr",
            ["NG_012337.1(NM_003002.2):c.(274G>T)"],
        ),
        (
            "NG_012337.1(NP_002993.1):p.Asp92Glu",
            [
                "NG_012337.1(NM_003002.2):c.(276C>G)",
                "NG_012337.1(NM_003002.2):c.(276C>A)",
            ],
        ),
        (
            "NG_012337.1(NP_002993.1):p.[Asp92Tyr;97Thr]",
            ["NG_012337.1(NM_003002.2):c.([274G>T;289G>A])"],
        ),
        (
            "NG_012337.1(NP_002993.1):p.[Asp92Glu;97Thr]",
            [
                "NG_012337.1(NM_003002.2):c.([276C>G;289G>A])",
                "NG_012337.1(NM_003002.2):c.([276C>A;289G>A])",
            ],
        ),
    ],
)
def test_protein(protein_description, back_translated_descriptions):
    assert set(back_translate(protein_description)) == set(back_translated_descriptions)
