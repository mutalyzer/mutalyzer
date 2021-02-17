import bisect
import copy

from mutalyzer_crossmapper import Coding, Genomic
from mutalyzer_mutator.util import reverse_complement

from normalizer.util import get_end, get_start

from .to_hgvs_coordinates import genomic_to_point, reverse_strand_shift


def _get_downstream_cds_position(selector_model, crossmap, offset):
    cds_end_coordinate = crossmap.genomic_to_coordinate(selector_model["cds"][0][1])
    cds_end_coding = crossmap.coordinate_to_coding(cds_end_coordinate)
    return genomic_to_point(crossmap.genomic_to_coordinate(cds_end_coding[0] + offset))


def _point_to_cds_coordinate(point, selector_model, crossmap, end_location=False):
    genomic_to_coordinate = Genomic().genomic_to_coordinate
    point_position = point["position"]
    if selector_model.get("inverted"):
        # TODO: This should be checked in more detail.
        point_position -= 1
        if point.get("shift"):
            point_position -= point["shift"]
    coding = crossmap.coordinate_to_coding(point_position, degenerate=True)
    if coding[2] == 0:
        if end_location and coding[1] > 0:
            return genomic_to_point(genomic_to_coordinate(coding[0] + 1))
        return genomic_to_point(genomic_to_coordinate(coding[0]))
    elif coding[2] == 1:
        return _get_downstream_cds_position(selector_model, crossmap, coding[0])
    elif coding[2] == -1:
        return genomic_to_point(0)
    # TODO: Should not reach this, but what if?


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
        if insertion.get("inverted"):
            inserted_value = reverse_complement(inserted_value)

    return {"source": "description", "sequence": inserted_value}


def variant_to_cds_coordinate(variant, sequences, selector_model, crossmap):
    new_variant = copy.deepcopy(variant)

    location = new_variant["location"]

    if location["type"] == "range":
        location["start"] = _point_to_cds_coordinate(
            location["start"], selector_model, crossmap
        )
        location["end"] = _point_to_cds_coordinate(
            location["end"], selector_model, crossmap, True
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


def location_in_same_intron(location, selector_model, coding_exons=True):
    exons = [e for l in selector_model["exon"] for e in l]
    cds = selector_model["cds"][0]
    if coding_exons:
        l_index = bisect.bisect_right(exons, cds[0])
        r_index = bisect.bisect_left(exons, cds[1])
        exons = [cds[0]] + exons[l_index:r_index] + [cds[1]]

    start_i = bisect.bisect_right(exons, get_start(location))
    end_i = bisect.bisect_left(exons, get_end(location))

    if start_i == end_i and start_i % 2 == 0:
        return True
    else:
        return False


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
        if (
            variant["type"] == "deletion_insertion"
            and variant.get("location")
            and not location_in_same_intron(variant["location"], selector_model)
        ):
            new_variants.append(
                variant_to_cds_coordinate(variant, sequences, selector_model, crossmap)
            )
    return new_variants
