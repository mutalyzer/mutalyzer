import copy

from mutalyzer_crossmapper import Coding, Genomic, NonCoding
from mutalyzer_mutator.util import reverse_complement

from ..description_model import (
    get_reference_id,
    yield_point_locations_for_main_reference,
    yield_ranges_main_reference,
)
from ..reference import (
    get_coordinate_system_from_selector_id,
    get_internal_selector_model,
)
from ..util import get_start, set_by_path
from .to_hgvs_indexing import to_hgvs_indexing
from .to_internal_coordinates import get_coordinate_system


def genomic_to_point(genomic):
    point = {"type": "point", "position": genomic}
    return point


def coding_to_point(coding):
    position, offset, section = coding[:3]
    point = {"type": "point", "position": position}

    if section == -1:
        point["outside_cds"] = "upstream"
        point["position"] = abs(point["position"])
    elif section == 1:
        point["outside_cds"] = "downstream"

    if offset != 0:
        point["offset"] = {"value": offset}
    return point


def noncoding_to_point(noncoding):
    position, offset, section = noncoding[:3]
    point = {"type": "point", "position": position}

    if offset != 0:
        point["offset"] = {"value": offset}
    return point


def point_to_hgvs(
    point, crossmap_function, point_function, degenerate=False, inverted=False
):
    if point.get("uncertain"):
        return {"type": "point", "uncertain": True}
    else:
        point_position = point["position"]
        if inverted and point.get("shift"):
            point_position -= point["shift"]
        if degenerate:
            new_point = point_function(crossmap_function(point_position, degenerate))
        else:
            new_point = point_function(crossmap_function(point_position))
    if point.get("shift"):
        new_point["shift"] = point["shift"]
    return new_point


def crossmap_to_hgvs_setup(coordinate_system, selector_model=None, degenerate=False):
    """
    Returns a crossmap instance able to convert from the internal system
    to the to hgvs system.
    """
    if coordinate_system in ["g", "m", "p", None]:
        crossmap = Genomic()
        return {
            "crossmap_function": crossmap.coordinate_to_genomic,
            "point_function": genomic_to_point,
        }
    elif coordinate_system == "c":
        crossmap = Coding(
            selector_model["exon"], selector_model["cds"][0], selector_model["inverted"]
        )
        return {
            "crossmap_function": crossmap.coordinate_to_coding,
            "point_function": coding_to_point,
            "degenerate": degenerate,
            "inverted": selector_model["inverted"],
        }
    elif coordinate_system == "n":
        crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
        return {
            "crossmap_function": crossmap.coordinate_to_noncoding,
            "point_function": noncoding_to_point,
            "inverted": selector_model["inverted"],
        }
    else:
        raise Exception("Unsupported coordinate system: {}.".format(coordinate_system))

    return crossmap


def initialize_hgvs_model(internal_model, coordinate_system=None, selector_id=None):
    model = copy.deepcopy(internal_model)
    if coordinate_system:
        model["coordinate_system"] = coordinate_system
    if model.get("reference"):
        if selector_id:
            model["reference"]["selector"] = {"id": selector_id}
        else:
            model["reference"] = {"id": model["reference"]["id"]}
    return model


def locations_to_hgvs_locations(internal_model, crossmap):
    hgvs_model = copy.deepcopy(internal_model)

    for point, path in yield_point_locations_for_main_reference(internal_model):
        set_by_path(hgvs_model, path, point_to_hgvs(point, **crossmap))

    return hgvs_model


def reverse_strand_shift(variants, seq):
    for variant in variants:
        if variant.get("inserted"):
            variant["inserted"].reverse()
            if (
                len(variant["inserted"]) == 1
                and variant["inserted"][0].get("sequence")
                and variant["location"]["start"].get("shift")
            ):
                # TODO: Check what to do when there is a compound insertion with locations included.
                start = get_start(variant)
                shift = variant["location"]["start"]["shift"]
                ins_seq = variant["inserted"][0]["sequence"]
                new_ins_seq = reverse_complement(
                    (seq[start - shift : start] + ins_seq)[: len(ins_seq)]
                )
                variant["inserted"][0]["sequence"] = new_ins_seq
            else:
                for inserted in variant["inserted"]:
                    if inserted.get("sequence"):
                        inserted["sequence"] = reverse_complement(inserted["sequence"])
        if variant.get("deleted"):
            variant["deleted"].reverse()
            for deleted in variant["deleted"]:
                if deleted.get("sequence"):
                    deleted["sequence"] = reverse_complement(deleted["sequence"])


def to_hgvs_locations(
    model,
    references,
    to_coordinate_system=None,
    to_selector_id=None,
    degenerate=False,
    selector_model=None,
):
    reference_id = get_reference_id(model)

    if to_selector_id and selector_model is None:
        selector_model = get_internal_selector_model(
            references[reference_id]["annotations"], to_selector_id, True
        )
    if to_selector_id is None and selector_model:
        to_selector_id = selector_model["id"]

    if not to_coordinate_system and selector_model:
        to_coordinate_system = get_coordinate_system_from_selector_id(
            references[reference_id], to_selector_id
        )

    hgvs_model = initialize_hgvs_model(model, to_coordinate_system, to_selector_id)

    to_coordinate_system = get_coordinate_system(hgvs_model, references)
    crossmap = crossmap_to_hgvs_setup(to_coordinate_system, selector_model, degenerate)

    if selector_model and selector_model.get("inverted"):
        reverse_strand_shift(
            hgvs_model["variants"], references["reference"]["sequence"]["seq"]
        )

    model_internal = to_hgvs_indexing(hgvs_model)

    for point, path in yield_point_locations_for_main_reference(model_internal):
        set_by_path(hgvs_model, path, point_to_hgvs(point, **crossmap))

    if selector_model and selector_model.get("inverted"):
        if hgvs_model.get("variants"):
            hgvs_model["variants"].reverse()
        for range_location, path in yield_ranges_main_reference(hgvs_model):
            range_location["start"], range_location["end"] = (
                range_location["end"],
                range_location["start"],
            )
    return hgvs_model
