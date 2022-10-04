import pytest

from mutalyzer.normalizer import normalize
from mutalyzer.reference import retrieve_reference, slice_to_selector

from .commons import code_in, patch_retriever
from .variants_set import TESTS_ALL


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test") and test.get(t_type):
            output.append((test["input"], test[t_type]))
    return output


@pytest.mark.parametrize(
    "input_description, normalized",
    get_tests([t for t in TESTS_ALL if "ensembl" in t["keywords"]], "normalized"),
)
def test_ensembl(input_description, normalized):
    d = normalize(input_description)
    assert d["normalized_description"] == normalized


@pytest.mark.parametrize(
    "id_ncbi, id_ensembl",
    [
        ("NM_003002.4", "ENST00000375549.8"),
        ("NM_024426.6", "ENST00000452863.10"),
        ("NM_003002.4", "ENST00000375549.8"),
    ],
)
def test_ensembl_mane_transcript(id_ncbi, id_ensembl):

    m_ncbi = retrieve_reference(id_ncbi)[0]
    m_ensembl = retrieve_reference(id_ensembl)[0]

    assert m_ncbi["sequence"]["seq"] == slice_to_selector(m_ensembl, id_ensembl, True)


@pytest.mark.parametrize(
    "id_ncbi, id_ensembl_gene, id_ensembl_transcript",
    [
        ("NM_003002.4", "ENSG00000204370.13", "ENST00000375549.8"),
        ("NM_024426.6", "ENSG00000184937.16", "ENST00000452863.10"),
    ],
)
def test_ensembl_mane_gene_transcript(id_ncbi, id_ensembl_gene, id_ensembl_transcript):
    m_ncbi = retrieve_reference(id_ncbi)[0]
    m_ensembl = retrieve_reference(id_ensembl_gene)[0]

    assert m_ncbi["sequence"]["seq"] == slice_to_selector(
        m_ensembl, id_ensembl_transcript, True
    )
