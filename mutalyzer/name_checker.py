from algebra import Variant
from algebra.extractor import extract_supremal, to_hgvs
from algebra.relations.supremal_based import find_supremal, spanning_variant
from algebra.variants import patch
from mutalyzer_hgvs_parser import to_model

from mutalyzer.util import get_end, get_inserted_sequence, get_start

from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .description import Description


def name_check_alt(description, only_variants=False, sequence=None):

    # TODO: reverse strand shift?

    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )

    d.to_delins()

    if d.errors:
        return d.output()

    if not only_variants and d.corrected_model["type"] == "description_protein":
        return {
            "errors": [
                {
                    "code": "ENOPROTEINSUPPORT",
                    "details": f"Protein descriptions not supported in this experimental service.",
                }
            ]
        }

    variants_algebra = []
    for variant in d.delins_model["variants"]:
        variants_algebra.append(
            Variant(
                get_start(variant),
                get_end(variant),
                get_inserted_sequence(variant, d.get_sequences()),
            )
        )

    reference = d.references["reference"]["sequence"]["seq"]
    observed = patch(reference, variants_algebra)
    supremal = find_supremal(
        reference, spanning_variant(reference, observed, variants_algebra)
    )
    canonical = list(extract_supremal(reference, supremal))
    algebra_hgvs = to_hgvs(canonical, reference)

    if only_variants:
        d.normalized_description = algebra_hgvs
        d.de_hgvs_model = {"variants": to_model(algebra_hgvs, "variants")}
        output = d.output()
    else:
        algebra_model = {
            "type": d.corrected_model["type"],
            "reference": {"id": d.corrected_model["reference"]["id"]},
            "coordinate_system": "g",
            "variants": to_model(algebra_hgvs, "variants"),
        }

        d.de_hgvs_internal_indexing_model = to_internal_indexing(
            to_internal_coordinates(algebra_model, d.get_sequences())
        )
        d.construct_de_hgvs_internal_indexing_model()
        d.construct_de_hgvs_coordinates_model()
        d.construct_normalized_description()
        d.construct_equivalent()

        output = d.output()
        output["algebra"] = algebra_hgvs

    output["supremal"] = {
        "hgvs": f"{d.corrected_model['reference']['id']}:g.{supremal.to_hgvs()}",
        "spdi": supremal.to_spdi(d.corrected_model["reference"]["id"])}

    return output


def name_check(description, only_variants=False, sequence=None):
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )

    d.normalize()

    return d.output()
