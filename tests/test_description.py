import pytest

from mutalyzer.description import Description

ASSEMBLY_TESTS = [
    (
        ["GRCh37", "hg19"],
        ["chr23", "23", "chrx", "x", "chrX", "X"],
        {"id": "NC_000023.10"},
    ),
    (
        ["GRCh37", "hg19"],
        ["chr24", "24", "chry", "y", "chrY", "Y"],
        {"id": "NC_000024.9"},
    ),
    (
        ["GRCh37", "hg19"],
        ["chr23", "23", "chrx", "x", "chrX", "X"],
        {"id": "NC_000023.10", "selector": {"id": "NM_004042.4"}},
        "NM_004042.4",
    ),
    (
        ["GRCh38", "hg38"],
        ["chr23", "23", "chrx", "x", "chrX", "X"],
        {"id": "NC_000023.11"},
    ),
    (
        ["GRCh38", "hg38"],
        ["chr24", "24", "chry", "y", "chrY", "Y"],
        {"id": "NC_000024.10"},
    ),
    (
        ["GRCh38", "hg38"],
        ["chr23", "23", "chrx", "x", "chrX", "X"],
        {"id": "NC_000023.11", "selector": {"id": "NM_004042.4"}},
        "NM_004042.4",
    ),
]


def get_tests():
    def _tests(assemblies, chromosomes, result, transcript=None):
        for assembly in assemblies:
            for chromosome in chromosomes:
                if transcript is None:
                    yield f"{assembly}({chromosome}):g.3250del", result
                else:
                    yield f"{assembly}({chromosome}({transcript})):g.3250del", result

    for t in ASSEMBLY_TESTS:
        yield from _tests(*t)


@pytest.mark.parametrize("description, reference_part", get_tests())
def test_assembly_chromosome_to_id(description, reference_part):
    d = Description(description)
    d.assembly_checks()
    assert d.corrected_model["reference"] == reference_part
