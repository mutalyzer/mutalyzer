import pytest

from normalizer.viewer import view_variants

from .commons import patch_retriever


@pytest.mark.parametrize(
    "description, output",
    [
        (
            "NM_003002.2:r.274g>u",
            [
                {
                    "description": "274g>u",
                    "seq_length": 1382,
                    "left": "ugaauccuugcucugcgaug",
                    "deleted": {"seq": "g"},
                    "inserted": {"seq": "u"},
                    "right": "acuauucccuggcugcagcc",
                    "start": 314,
                }
            ],
        )
    ],
)
def test_view_variants(description, output):
    assert view_variants(description) == output
