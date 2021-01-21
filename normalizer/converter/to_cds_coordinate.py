import copy
import json

from mutalyzer_crossmapper import Coding, Genomic

from normalizer.util import get_end, get_start

from .to_hgvs import genomic_to_point
from .to_hgvs_coordinates import reverse_strand_shift


def _get_last_exon_cds_coding(selector_model, crossmap):
    last_exon_coordinate = crossmap.genomic_to_coordinate(
        selector_model["exon"][-1][-1]
    )
    last_exon_coding = crossmap.coordinate_to_coding(last_exon_coordinate)
    return genomic_to_point(crossmap.genomic_to_coordinate(last_exon_coding[0]))


def _point_to_cds_coordinate(point, selector_model, crossmap):
    genomic_to_coordinate = Genomic().genomic_to_coordinate
    point_position = point["position"]
    if selector_model.get("inverted"):
        # TODO: This should be checked in more detail.
        point_position -= 1
        if point.get("shift"):
            point_position -= point["shift"]
    coding = crossmap.coordinate_to_coding(point_position)
    if coding[2] == 0:
        return genomic_to_point(genomic_to_coordinate(coding[0]))
    elif coding[2] == -1:
        return genomic_to_point(genomic_to_coordinate(coding[0]))
    elif coding[2] == 1:
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
    new_variant["location"] = location
    return new_variant


def reverse_start_end(variants):
    for variant in variants:
        if variant.get("location") and variant["location"]["type"] == "range":
            loc = variant["location"]
            loc["start"], loc["end"] = loc["end"], loc["start"]


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
    shifted_variants = copy.deepcopy(variants)
    if selector_model.get("inverted"):
        reverse_strand_shift(shifted_variants, sequences["reference"])
        reverse_start_end(shifted_variants)

    new_variants = []
    for variant in shifted_variants:
        if variant["type"] == "deletion_insertion":
            new_variants.append(
                variant_to_cds_coordinate(variant, sequences, selector_model, crossmap)
            )
    return new_variants
