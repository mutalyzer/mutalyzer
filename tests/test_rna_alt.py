import pytest

from mutalyzer.rna import dna_to_rna, rna_to_dna, get_position_type
from mutalyzer.description_model import model_to_string

from .commons import code_in, monkey_patches
from .variants_set import TESTS_ALL

EXONS = [(5026, 5113), (6010, 6127), (7020, 7165), (12958, 13948)]


@pytest.mark.parametrize(
    "position, exons, len_ss, position_type",
    [
        # in intron
        (7015, EXONS, 4, (4, 0)),
        (7016, EXONS, 4, (4, -4)),
        (7017, EXONS, 4, (4, -3)),
        (7018, EXONS, 4, (4, -2)),
        (7019, EXONS, 4, (4, -1)),
        # in exon
        (7020, EXONS, 2, (5, 1)),
        (7021, EXONS, 2, (5, 2)),
        (7022, EXONS, 2, (5, 0)),
        (7162, EXONS, 2, (5, 0)),
        (7163, EXONS, 2, (5, -2)),
        (7164, EXONS, 2, (5, -1)),
        # in intron
        (7165, EXONS, 4, (6, 1)),
        (7166, EXONS, 4, (6, 2)),
        (7167, EXONS, 4, (6, 3)),
        (7168, EXONS, 4, (6, 4)),
        (7169, EXONS, 4, (6, 0)),
    ],
)
def test_get_position_type(position, exons, len_ss, position_type):
    assert get_position_type(position, exons, len_ss) == position_type


def get_tests():
    output = []
    for test in TESTS_ALL:
        if test.get("to_test"):
            if test.get("rna_description_alt") is not None:
                if test.get("rna_description_alt") is False:
                    output.append((test["input"], None))
                else:
                    output.append((test["input"], test.get("rna_description_alt")))
            elif test.get("rna_description"):
                output.append((test["input"], test["rna_description"]))
    return output


@pytest.mark.parametrize(
    "input_description, rna_expected",
    get_tests(),
)
def test_dna_to_rna_variants_set(input_description, rna_expected):
    rna = dna_to_rna(input_description)
    assert rna.get("description") == rna_expected


@pytest.mark.parametrize(
    "input_description, rna_expected",
    [
        # First exon start: c.-35
        ("NG_012337.3(NM_003002.4):c.-41_-40insA", "NG_012337.3(NM_003002.4):r.(=)"),
        ("NG_012337.3(NM_003002.4):c.-40_-39insA", "NG_012337.3(NM_003002.4):r.(=)"),
        ("NG_012337.3(NM_003002.4):c.-34_-33insA", "NG_012337.3(NM_003002.4):r.(-34_-33insa)"),

        # First exon end: c.52
        ("NG_012337.3(NM_003002.4):c.50_51insC", "NG_012337.3(NM_003002.4):r.(50_51insc)"),
        ("NG_012337.3(NM_003002.4):c.52+4_52+5insC", "NG_012337.3(NM_003002.4):r.(=)"),

        # Exon 2 start: c.53
        ("NG_012337.3(NM_003002.4):c.53-5_53-4insA", "NG_012337.3(NM_003002.4):r.(=)"),
        ("NG_012337.3(NM_003002.4):c.54_55insA", "NG_012337.3(NM_003002.4):r.(54_55insa)"),

        # Exon 2 end: c.169
        ("NG_012337.3(NM_003002.4):c.167_168insC", "NG_012337.3(NM_003002.4):r.(167_168insc)"),
        ("NG_012337.3(NM_003002.4):c.169+4_169+5insC", "NG_012337.3(NM_003002.4):r.(=)"),

        # Last exon start: c.315
        ("NG_012337.3(NM_003002.4):c.315-5_315-4insC", "NG_012337.3(NM_003002.4):r.(=)"),
        ("NG_012337.3(NM_003002.4):c.316_317insC", "NG_012337.3(NM_003002.4):r.(316_317insc)"),

        # Last exon end: c.*824
        ("NG_012337.3(NM_003002.4):c.*821_*822insC", "NG_012337.3(NM_003002.4):r.(*821_*822insc)"),
        ("NG_012337.3(NM_003002.4):c.*822_*823insC", "NG_012337.3(NM_003002.4):r.(*822_*823insc)"),
        ("NG_012337.3(NM_003002.4):c.*828_*829insC", "NG_012337.3(NM_003002.4):r.(=)"),

        # Other
        ("NG_012337.1(NM_012459.2):c.271del", "NG_012337.1(NM_012459.2):r.(269_271c[2])"),
        ("NM_012459.2:c.271del", "NM_012459.2:r.(269_271c[2])"),
        ("NG_012337.1(NM_012459.2):c.271C>A", "NG_012337.1(NM_012459.2):r.(271c>a)"),
        ("NG_012337.3(NM_003002.4):c.274G>T", "NG_012337.3(NM_003002.4):r.(274g>u)"),
        ("NM_003002.4:c.274G>T", "NM_003002.4:r.(274g>u)"),
        ("NG_012337.1(NM_003002.2):c.166C>A", "NG_012337.1(NM_003002.2):r.(166c>a)"),
        ("NR_038420.1:n.206_210del", "NR_038420.1:r.(206_210del)"),
        ("NG_007485.1(NR_024274.1):n.211delinsGG", "NG_007485.1(NR_024274.1):r.(211_215g[6])"),
        ("NM_003002.4:c.169_170insAAA", "NM_003002.4:r.(169_170insa[3])"),
    ],
)
def test_dna_to_rna_new(input_description, rna_expected):
    rna = dna_to_rna(input_description)
    assert rna.get("description") == rna_expected


@pytest.mark.parametrize(
    "description",
    [
        # First exon start: c.-35
        "NG_012337.3(NM_003002.4):c.-39_-38insA",
        "NG_012337.3(NM_003002.4):c.-38_-37insA",
        "NG_012337.3(NM_003002.4):c.-37_-36insA",
        "NG_012337.3(NM_003002.4):c.-36_-35insA",
        "NG_012337.3(NM_003002.4):c.-35_-34insA",

        # First exon end: c.52
        "NG_012337.3(NM_003002.4):c.51_52insC",
        "NG_012337.3(NM_003002.4):c.52_52+1insC",
        "NG_012337.3(NM_003002.4):c.52+1_52+2insC",
        "NG_012337.3(NM_003002.4):c.52+2_52+3insC",
        "NG_012337.3(NM_003002.4):c.52+3_52+4insC",

        # Exon 2 start: c.53
        "NG_012337.3(NM_003002.4):c.53-4_53-3insA",
        "NG_012337.3(NM_003002.4):c.53-3_53-2insA",
        "NG_012337.3(NM_003002.4):c.53-2_53-1insA",
        "NG_012337.3(NM_003002.4):c.53-1_53insA",
        "NG_012337.3(NM_003002.4):c.53_54insA",

        # Exon 2 end: c.169
        "NG_012337.3(NM_003002.4):c.168_169insC",
        "NG_012337.3(NM_003002.4):c.169_169+1insC",
        "NG_012337.3(NM_003002.4):c.169+1_169+2insC",
        "NG_012337.3(NM_003002.4):c.169+2_169+3insC",
        "NG_012337.3(NM_003002.4):c.169+3_169+4insC",

        # Last exon start: c.315
        "NG_012337.3(NM_003002.4):c.315-4_315-3insC",
        "NG_012337.3(NM_003002.4):c.315-3_315-2insC",
        "NG_012337.3(NM_003002.4):c.315-2_315-1insC",
        "NG_012337.3(NM_003002.4):c.315-1_315insC",
        "NG_012337.3(NM_003002.4):c.315_316insC",

        # Last exon end: c.*824
        "NG_012337.3(NM_003002.4):c.*823_*824insC",
        "NG_012337.3(NM_003002.4):c.*823del",
        "NG_012337.3(NM_003002.4):c.*824del",
        "NG_012337.3(NM_003002.4):c.*824_*825insC",
        "NG_012337.3(NM_003002.4):c.*825_*826insC",
        "NG_012337.3(NM_003002.4):c.*826_*827insC",
        "NG_012337.3(NM_003002.4):c.*827_*828insC",

        # Other
        "NG_012337.1(NM_003002.2):c.168del",
        "NG_012337.1(NM_003002.2):c.169del",
        "NG_012337.1(NM_003002.2):c.170del",
        "NG_012337.3(NM_003002.4):c.169_169+1insAA",
        "NG_012337.3(NM_003002.4):c.169_170insAA",
        "NG_007485.1(NR_024274.1):n.616_616+1insTTTTTT",
        "NG_012337.3(NM_003002.4):c.53-20_169+10del",
        "NG_008835.1(NM_022153.2):c.677-20_704+62del",
        # "NG_008835.1(NM_022153.2):c.83-20_511+62del",
        "NG_008835.1(NM_022153.2):c.512-20_568+62del",
    ],
)
def test_dna_to_rna_new_errors(description):
    assert dna_to_rna(description).get("errors") is not None


@pytest.mark.parametrize(
    "input_description, dna_expected",
    [
        ("NG_012337.1(NM_012459.2):r.(269_271c[2])", "NG_012337.1(NM_012459.2):c.(269_271C[2])"),
        ("NM_012459.2:r.(269_271c[2])", "NM_012459.2:c.(269_271C[2])"),
        ("NG_012337.1(NM_012459.2):r.(271c>a)", "NG_012337.1(NM_012459.2):c.(271C>A)"),
        ("NG_012337.3(NM_003002.4):r.(274g>u)", "NG_012337.3(NM_003002.4):c.(274G>T)"),
        ("NM_003002.4:r.(274g>u)", "NM_003002.4:c.(274G>T)"),
        ("NG_012337.1(NM_003002.2):r.(166c>a)", "NG_012337.1(NM_003002.2):c.(166C>A)"),
        ("NR_038420.1:r.(206_210del)", "NR_038420.1:n.(206_210del)"),
        ("NG_007485.1(NR_024274.1):r.(211_215g[6])", "NG_007485.1(NR_024274.1):n.(211_215G[6])"),
        ("NM_003002.4:r.(169_170insa[3])", "NM_003002.4:c.(169_170insA[3])"),
    ],
)
def test_rna_to_dna_new(input_description, dna_expected):
    dna = rna_to_dna(input_description)
    assert model_to_string(dna) == dna_expected


@pytest.mark.parametrize(
    "description",
    [
        "NG_012337.3(NM_003002.4):r.169_170insAAA",
        "NG_007485.1(NR_024274.1):r.616_617insTTTTTT",
    ],
)
def test_rna_to_dna_new_errors(description):
    assert rna_to_dna(description).get("errors") is not None
