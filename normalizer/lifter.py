import extractor

from .converter import de_to_hgvs
from .converter.to_hgvs_coordinates import (
    crossmap_to_hgvs_setup,
    locations_to_hgvs_locations,
    to_hgvs_locations,
)
from .converter.to_hgvs_indexing import variants_to_internal_indexing
from .description import Description
from .description_model import model_to_string, variants_to_description
from .reference import get_reference_model, get_reference_mol_type


def lift(description, reference_id):
    d = Description(description)
    d.normalize()
    if d.references and d.references.get("observed"):
        r_model = get_reference_model(reference_id)
        ref_seq = r_model["sequence"]["seq"]
        obs_seq = d.references["observed"]["sequence"]["seq"]
        de_hgvs_variants = de_to_hgvs(
            extractor.describe_dna(ref_seq, obs_seq),
            {"reference": ref_seq, "observed": obs_seq},
        )
        if get_reference_mol_type(r_model) in ["mRNA"]:
            return model_to_string(
                to_hgvs_locations(
                    model={
                        "reference": {"id": reference_id},
                        "variants": de_hgvs_variants,
                    },
                    references={"reference": r_model, reference_id: r_model},
                    to_coordinate_system="c",
                    to_selector_id=reference_id,
                    degenerate=True,
                )
            )
        else:
            hgvs_indexing_variants = variants_to_internal_indexing(de_hgvs_variants)
            crossmap = crossmap_to_hgvs_setup("g")
            hgvs_variants = locations_to_hgvs_locations(
                {"variants": hgvs_indexing_variants}, crossmap
            )
            return model_to_string(
                {
                    "reference": {"id": reference_id},
                    "coordinate_system": "g",
                    "variants": hgvs_variants["variants"],
                }
            )
