import copy
import json

from mutalyzer_crossmapper import Genomic, Coding

from .to_hgvs import genomic_to_point
from normalizer.util import get_start, get_end


def _get_last_exon_cds_coding(selector_model, crossmap):
    last_exon_coordinate = crossmap.genomic_to_coordinate(
        selector_model["exon"][-1][-1]
    )
    last_exon_coding = crossmap.coordinate_to_coding(last_exon_coordinate)
    return genomic_to_point(crossmap.genomic_to_coordinate(last_exon_coding[0]))


def _point_to_cds_coordinate(point, selector_model, crossmap):
    genomic_to_coordinate = Genomic().genomic_to_coordinate
    coding = crossmap.coordinate_to_coding(point["position"])
    if coding[2] == 0:
        return genomic_to_point(0)
    elif coding[2] == 1:
        return genomic_to_point(genomic_to_coordinate(coding[0]))
    elif coding[2] == 2:
        return _get_last_exon_cds_coding(selector_model, crossmap)


def get_inserted_sequence(insertion, sequences):
    return sequences[insertion["source"]][
        get_start(insertion["location"]) : get_end(insertion["location"])
    ]


def merge_inserted_to_string(inserted, sequences):
    inserted_value = ""
    for insertion in inserted:
        if insertion.get("sequence"):
            inserted_value += insertion.get("sequence")
        else:
            inserted_value += get_inserted_sequence(insertion, sequences)
    return {"source": "description", "sequence": inserted_value}


def variant_to_cds_coordinate(variant, sequences, selector_model, crossmap):
    new_variant = copy.deepcopy(variant)

    location = new_variant["location"]

    if location["type"] == "range":
        location["start"] = _point_to_cds_coordinate(
            location["start"], selector_model, crossmap
        )
        location["end"] = _point_to_cds_coordinate(
            location["end"], selector_model, crossmap
        )
    else:
        location = _point_to_cds_coordinate(location, selector_model, crossmap)
    if new_variant.get("inserted"):
        new_variant["inserted"] = [
            merge_inserted_to_string(new_variant["inserted"], sequences)
        ]
    return new_variant


def to_cds_coordinate(variants, sequences, selector_model):
    """
    Converts the locations to cds equivalent.

    :param variants: Variants with locations in the coordinate system.
    :param selector_model:
    :param crossmap:
    """
    crossmap = Coding(
        selector_model["exon"], selector_model["cds"][0], selector_model["inverted"]
    )
    new_variants = []
    for variant in variants:
        if variant["type"] == "deletion_insertion":
            new_variants.append(
                variant_to_cds_coordinate(variant, sequences, selector_model, crossmap)
            )
    return new_variants
