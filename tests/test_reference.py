import pytest

from mutalyzer.reference import (
    get_internal_selector_model,
    get_reference_model_segmented,
    get_selector_feature,
)

from .commons import patch_retriever


def test_get_reference_model_segmented_record_no_siblings_ancestors_no_descendants():
    feature_model = {
        "id": "NG_012337.1",
        "type": "record",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": 15948},
        },
        "qualifiers": {
            "mol_type": "genomic DNA",
            "chromosome": "11",
            "dbxref": "taxon:9606",
            "map": "11q23",
            "name": "11",
            "genome": "genomic",
        },
    }
    assert feature_model == get_reference_model_segmented(
        "NG_012337.1", "NG_012337.1", False, True, False
    )


def test_get_reference_model_segmented_transcript_no_siblings_ancestors_descendants():
    feature_model = {
        "id": "NG_012337.1",
        "type": "record",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": 15948},
        },
        "qualifiers": {
            "mol_type": "genomic DNA",
            "chromosome": "11",
            "dbxref": "taxon:9606",
            "map": "11q23",
            "name": "11",
            "genome": "genomic",
        },
        "features": [
            {
                "id": "C11orf57",
                "type": "gene",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 0},
                    "end": {"type": "point", "position": 3304},
                    "strand": 1,
                },
                "qualifiers": {"name": "C11orf57", "HGNC": "25569"},
                "features": [
                    {
                        "id": "NM_018195.3",
                        "type": "mRNA",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 135},
                            "end": {"type": "point", "position": 3304},
                            "strand": 1,
                        },
                        "features": [
                            {
                                "id": "exon-NM_018195.3-1",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 135},
                                    "end": {"type": "point", "position": 189},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "exon-NM_018195.3-2",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 618},
                                    "end": {"type": "point", "position": 3304},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "NP_060665.3",
                                "type": "CDS",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 135},
                                    "end": {"type": "point", "position": 1126},
                                    "strand": 1,
                                },
                                "qualifiers": {
                                    "exception": "annotated by transcript or proteomic data"
                                },
                            },
                        ],
                    },
                ],
            },
        ],
    }
    assert feature_model == get_reference_model_segmented("NG_012337.1", "NM_018195.3")


def test_get_reference_model_segmented_transcript_no_siblings_ancestors_no_descendants():
    feature_model = {
        "id": "NG_012337.1",
        "type": "record",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": 15948},
        },
        "qualifiers": {
            "mol_type": "genomic DNA",
            "chromosome": "11",
            "dbxref": "taxon:9606",
            "map": "11q23",
            "name": "11",
            "genome": "genomic",
        },
        "features": [
            {
                "id": "C11orf57",
                "type": "gene",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 0},
                    "end": {"type": "point", "position": 3304},
                    "strand": 1,
                },
                "qualifiers": {"name": "C11orf57", "HGNC": "25569"},
                "features": [
                    {
                        "id": "NM_018195.3",
                        "type": "mRNA",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 135},
                            "end": {"type": "point", "position": 3304},
                            "strand": 1,
                        },
                    },
                ],
            },
        ],
    }
    assert feature_model == get_reference_model_segmented(
        "NG_012337.1", "NM_018195.3", False, True, False
    )


def test_get_reference_model_segmented_transcript_no_siblings_no_ancestors_descendants():
    feature_model = {
        "id": "NM_018195.3",
        "type": "mRNA",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 135},
            "end": {"type": "point", "position": 3304},
            "strand": 1,
        },
        "features": [
            {
                "id": "exon-NM_018195.3-1",
                "type": "exon",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 135},
                    "end": {"type": "point", "position": 189},
                    "strand": 1,
                },
            },
            {
                "id": "exon-NM_018195.3-2",
                "type": "exon",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 618},
                    "end": {"type": "point", "position": 3304},
                    "strand": 1,
                },
            },
            {
                "id": "NP_060665.3",
                "type": "CDS",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 135},
                    "end": {"type": "point", "position": 1126},
                    "strand": 1,
                },
                "qualifiers": {
                    "exception": "annotated by transcript or proteomic data"
                },
            },
        ],
    }
    assert feature_model == get_reference_model_segmented(
        "NG_012337.1", "NM_018195.3", False, False, True
    )


def test_get_reference_model_segmented_transcript_no_siblings_no_descendants_no_ancestors():
    feature_model = {
        "id": "NM_018195.3",
        "type": "mRNA",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 135},
            "end": {"type": "point", "position": 3304},
            "strand": 1,
        },
    }
    assert feature_model == get_reference_model_segmented(
        "NG_012337.1", "NM_018195.3", False, False, False
    )


def test_get_reference_model_segmented_transcript_include_siblings():
    feature_model = {
        "id": "NG_012337.1",
        "type": "record",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": 15948},
        },
        "qualifiers": {
            "mol_type": "genomic DNA",
            "chromosome": "11",
            "dbxref": "taxon:9606",
            "map": "11q23",
            "name": "11",
            "genome": "genomic",
        },
        "features": [
            {
                "id": "SDHD",
                "type": "gene",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 5000},
                    "end": {"type": "point", "position": 13948},
                    "strand": 1,
                },
                "qualifiers": {
                    "name": "SDHD",
                    "synonym": ["CBT1", "PGL", "PGL1", "SDH4"],
                    "HGNC": "10683",
                },
                "features": [
                    {
                        "id": "id-SDHD-1",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 5000},
                            "end": {"type": "point", "position": 5113},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "id-SDHD-2",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 6010},
                            "end": {"type": "point", "position": 6127},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "id-SDHD-3",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 7020},
                            "end": {"type": "point", "position": 7165},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "id-SDHD-4",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 12958},
                            "end": {"type": "point", "position": 13948},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "NM_003002.2",
                        "type": "mRNA",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 5000},
                            "end": {"type": "point", "position": 13948},
                            "strand": 1,
                        },
                        "features": [
                            {
                                "id": "exon-NM_003002.2-1",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 5000},
                                    "end": {"type": "point", "position": 5113},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "exon-NM_003002.2-2",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 6010},
                                    "end": {"type": "point", "position": 6127},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "exon-NM_003002.2-3",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 7020},
                                    "end": {"type": "point", "position": 7165},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "exon-NM_003002.2-4",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 12958},
                                    "end": {"type": "point", "position": 13948},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "NP_002993.1",
                                "type": "CDS",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 5061},
                                    "end": {"type": "point", "position": 13124},
                                    "strand": 1,
                                },
                            },
                        ],
                    },
                ],
            }
        ],
    }

    assert feature_model == get_reference_model_segmented(
        "NG_012337.1", "NM_003002.2", True
    )


def test_get_reference_model_segmented_gene_no_siblings():
    feature_model = {
        "id": "NG_012337.1",
        "type": "record",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": 0},
            "end": {"type": "point", "position": 15948},
        },
        "qualifiers": {
            "mol_type": "genomic DNA",
            "chromosome": "11",
            "dbxref": "taxon:9606",
            "map": "11q23",
            "name": "11",
            "genome": "genomic",
        },
        "features": [
            {
                "id": "SDHD",
                "type": "gene",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 5000},
                    "end": {"type": "point", "position": 13948},
                    "strand": 1,
                },
                "qualifiers": {
                    "name": "SDHD",
                    "synonym": ["CBT1", "PGL", "PGL1", "SDH4"],
                    "HGNC": "10683",
                },
                "features": [
                    {
                        "id": "id-SDHD-1",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 5000},
                            "end": {"type": "point", "position": 5113},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "id-SDHD-2",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 6010},
                            "end": {"type": "point", "position": 6127},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "id-SDHD-3",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 7020},
                            "end": {"type": "point", "position": 7165},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "id-SDHD-4",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 12958},
                            "end": {"type": "point", "position": 13948},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "NM_003002.2",
                        "type": "mRNA",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 5000},
                            "end": {"type": "point", "position": 13948},
                            "strand": 1,
                        },
                        "features": [
                            {
                                "id": "exon-NM_003002.2-1",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 5000},
                                    "end": {"type": "point", "position": 5113},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "exon-NM_003002.2-2",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 6010},
                                    "end": {"type": "point", "position": 6127},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "exon-NM_003002.2-3",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 7020},
                                    "end": {"type": "point", "position": 7165},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "exon-NM_003002.2-4",
                                "type": "exon",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 12958},
                                    "end": {"type": "point", "position": 13948},
                                    "strand": 1,
                                },
                            },
                            {
                                "id": "NP_002993.1",
                                "type": "CDS",
                                "location": {
                                    "type": "range",
                                    "start": {"type": "point", "position": 5061},
                                    "end": {"type": "point", "position": 13124},
                                    "strand": 1,
                                },
                            },
                        ],
                    },
                ],
            }
        ],
    }

    assert feature_model == get_reference_model_segmented("NG_012337.1", "SDHD", False)


@pytest.mark.parametrize(
    "reference_id, feature_id, feature_model",
    [
        (
            "NG_012337.1",
            "NM_003002.2",
            {
                "id": "NM_003002.2",
                "type": "mRNA",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 5000},
                    "end": {"type": "point", "position": 13948},
                    "strand": 1,
                },
                "features": [
                    {
                        "id": "exon-NM_003002.2-1",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 5000},
                            "end": {"type": "point", "position": 5113},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "exon-NM_003002.2-2",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 6010},
                            "end": {"type": "point", "position": 6127},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "exon-NM_003002.2-3",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 7020},
                            "end": {"type": "point", "position": 7165},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "exon-NM_003002.2-4",
                        "type": "exon",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 12958},
                            "end": {"type": "point", "position": 13948},
                            "strand": 1,
                        },
                    },
                    {
                        "id": "NP_002993.1",
                        "type": "CDS",
                        "location": {
                            "type": "range",
                            "start": {"type": "point", "position": 5061},
                            "end": {"type": "point", "position": 13124},
                            "strand": 1,
                        },
                    },
                ],
            },
        ),
        (
            "NG_012337.1",
            "NP_002993.1",
            {
                "id": "NP_002993.1",
                "type": "CDS",
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 5061},
                    "end": {"type": "point", "position": 13124},
                    "strand": 1,
                },
            },
        ),
    ],
)
def test_get_feature(reference_id, feature_id, feature_model):
    reference_model = get_reference_model_segmented(reference_id)["annotations"]
    assert get_selector_feature(reference_model, feature_id) == feature_model


@pytest.mark.parametrize(
    "reference_id, feature_id, feature_model",
    [
        (
            "NG_012337.1",
            "NM_003002.2",
            {
                "id": "NM_003002.2",
                "cds_id": "NP_002993.1",
                "type": "mRNA",
                "inverted": False,
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 5000},
                    "end": {"type": "point", "position": 13948},
                    "strand": 1,
                },
                "exon": [(5000, 5113), (6010, 6127), (7020, 7165), (12958, 13948)],
                "cds": [(5061, 13124)],
            },
        ),
        (
            "NG_012337.1",
            "NP_002993.1",
            {
                "id": "NP_002993.1",
                "mrna_id": "NM_003002.2",
                "type": "CDS",
                "inverted": False,
                "location": {
                    "type": "range",
                    "start": {"type": "point", "position": 5061},
                    "end": {"type": "point", "position": 13124},
                    "strand": 1,
                },
                "exon": [(5000, 5113), (6010, 6127), (7020, 7165), (12958, 13948)],
                "cds": [(5061, 13124)],
            },
        ),
    ],
)
def test_get_internal_selector_model(reference_id, feature_id, feature_model):
    reference_model = get_reference_model_segmented(reference_id)["annotations"]
    assert get_internal_selector_model(reference_model, feature_id) == feature_model
