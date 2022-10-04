import pytest

from mutalyzer.normalizer import normalize

from .commons import code_in, patch_retriever
from .variants_set import TESTS_ALL


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(t_type):
            output.append((test["input"], test[t_type]))
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
        assert d["equivalent_descriptions"]["g"][0] == genomic


@pytest.mark.parametrize("input_description, genomic", get_tests(TESTS_ALL, "genomic"))
def test_genomic(input_description, genomic):
    d = normalize(input_description)
    if d["equivalent_descriptions"].get("g"):
        assert d["equivalent_descriptions"]["g"][0] == genomic


@pytest.mark.parametrize(
    "input_description, coding", get_tests(TESTS_ALL, "coding_protein_descriptions")
)
def test_coding(input_description, coding):
    d = normalize(input_description)
    coding = [c[0] for c in coding]
    if d["equivalent_descriptions"].get("c"):
        name_check_coding = [c[0] for c in d["equivalent_descriptions"]["c"]]
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
    "input_description, coding_protein_descriptions",
    get_tests(TESTS_ALL, "coding_protein_descriptions"),
)
def test_protein_equivalent(input_description, coding_protein_descriptions):
    normalized_output = normalize(input_description)
    normalizer_descriptions = set(normalized_output["equivalent_descriptions"]["c"])

    assert coding_protein_descriptions.issubset(normalizer_descriptions)


@pytest.mark.parametrize(
    "input_description, rna_description",
    get_tests(TESTS_ALL, "rna_description"),
)
def test_rna(input_description, rna_description):

    normalized_output = normalize(input_description)
    normalizer_protein = normalized_output["rna"]["description"]

    assert normalizer_protein == rna_description


@pytest.mark.parametrize("input_description, codes", get_tests(TESTS_ALL, "errors"))
def test_errors(input_description, codes):
    assert codes == [error["code"] for error in normalize(input_description)["errors"]]


@pytest.mark.parametrize("input_description, codes", get_tests(TESTS_ALL, "infos"))
def test_infos(input_description, codes):
    assert codes == [info["code"] for info in normalize(input_description)["infos"]]


@pytest.mark.parametrize(
    "description, sequence, normalized",
    [
        ("2del", "AAAT", "3del"),
        ("[2del]", "AAAT", "3del"),
        ("[1del;2del]", "AAAT", "2_3del"),
        ("1_2insNG_012337.1:g.100", "AAAT", "1_2insT"),
        ("[9dup;14_15insCCTCT]", "CTCTCTCTCTCTCTTG", "10delinsCTCTCTC"),
    ],
)
def test_only_variants(description, sequence, normalized):
    assert (
            normalize(description, True, sequence)["normalized_description"] == normalized
    )
