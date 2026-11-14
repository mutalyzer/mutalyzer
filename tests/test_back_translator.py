import pytest
from mutalyzer_retriever.request import RequestErrors

from mutalyzer.back_translator import back_translate
from mutalyzer.normalizer import normalize

from .commons import code_in, monkey_patches


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
        (
            "NG_012337.1(NM_003002.2):p.[Asp92Glu;97Thr]",
            [
                "NG_012337.1(NM_003002.2):c.([276C>G;289G>A])",
                "NG_012337.1(NM_003002.2):c.([276C>A;289G>A])",
            ],
        ),
        # Gene with one CDS and no transcript: the gene is the c. selector.
        ("L78833.1(rpL21):p.Met1Val", ["L78833.1(rpL21):c.(1A>G)"]),
        # Bare CDS: no transcript, so no c. description.
        ("KJ082064.1(AHN50432.1):p.Arg11Gln", []),
    ],
)
def test_protein(protein_description, back_translated_descriptions):
    assert set(back_translate(protein_description)) == set(back_translated_descriptions)


def _lookup_fails(cds_id):
    raise RequestErrors([ConnectionError("NCBI datasets unavailable")])


@pytest.mark.parametrize(
    "description, get_cds_to_mrna",
    [
        # The related mRNA is found but cannot be retrieved.
        ("NP_002993.1:p.Asp92Tyr", lambda cds_id: ["NM_999999.9"]),
        # The related mRNA lookup fails.
        ("NP_002993.1:p.Asp92Tyr", _lookup_fails),
        # No related mRNA (the lookup is patched to None in commons).
        ("NP_002993.1:p.Asp92Tyr", None),
        # Bare CDS: no transcript in the annotation.
        ("KJ082064.1(AHN50432.1):p.Arg11Gln", None),
    ],
)
def test_protein_no_back_translation_info(description, get_cds_to_mrna, monkeypatch):
    if get_cds_to_mrna:
        monkeypatch.setattr("mutalyzer.description.get_cds_to_mrna", get_cds_to_mrna)
    output = normalize(description)
    assert "errors" not in output
    assert "back_translated_descriptions" not in output
    assert code_in("IBACKTRANSLATION", output["infos"])


def test_protein_selector_needs_no_lookup(monkeypatch):
    # The transcript is in the reference annotation: no lookup is needed.
    monkeypatch.setattr("mutalyzer.description.get_cds_to_mrna", _lookup_fails)
    output = normalize("NG_012337.1(NP_002993.1):p.Asp92Tyr")
    assert output["back_translated_descriptions"] == ["NG_012337.1(NM_003002.2):c.(274G>T)"]
    assert not code_in("IBACKTRANSLATION", output.get("infos", []))
