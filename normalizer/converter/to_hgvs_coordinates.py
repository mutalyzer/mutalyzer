import copy

from mutalyzer_crossmapper import Coding, Genomic, NonCoding

from ..description_model import (
    get_reference_id,
    yield_point_locations_for_main_reference,
)
from ..reference import get_coordinate_system_from_selector_id, get_selector_model
from ..util import set_by_path


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


def point_to_hgvs(point, crossmap_function, point_function, degenerate=False):
    if point.get("uncertain"):
        return {"type": "point", "uncertain": True}
    else:
        if degenerate:
            return point_function(crossmap_function(point["position"], degenerate))
        else:
            return point_function(crossmap_function(point["position"]))


def crossmap_to_hgvs_setup(coordinate_system, selector_model=None, degenerate=False):
    """
    Returns a crossmap instance able to convert from the internal system
    to the to hgvs system.
    """
    if coordinate_system == "g":
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
        }
    elif coordinate_system == "n":
        crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
        return {
            "crossmap_function": crossmap.coordinate_to_noncoding,
            "point_function": noncoding_to_point,
            "degenerate": degenerate,
        }
    else:
        raise Exception("Unsupported coordinate system: {}.".format(coordinate_system))

    return crossmap


def initialize_hgvs_model(internal_model, coordinate_system=None, selector_id=None):
    model = copy.deepcopy(internal_model)
    if coordinate_system:
        model["coordinate_system"] = coordinate_system
    if selector_id:
        model["reference"]["selector"] = {"id": selector_id}
    return model


def locations_to_hgvs_locations(internal_model, crossmap):
    hgvs_model = copy.deepcopy(internal_model)

    for point, path in yield_point_locations_for_main_reference(internal_model):
        set_by_path(hgvs_model, path, point_to_hgvs(point, **crossmap))

    return hgvs_model


def to_hgvs_locations(
    internal_model,
    references,
    to_coordinate_system=None,
    to_selector_id=None,
    degenerate=False,
):
    reference_id = get_reference_id(internal_model)

    selector_model = (
        get_selector_model(references[reference_id]["annotations"], to_selector_id)
        if to_selector_id
        else None
    )

    if to_coordinate_system is None and selector_model:
        to_coordinate_system = get_coordinate_system_from_selector_id(
            references[reference_id], to_selector_id
        )

    hgvs_model = initialize_hgvs_model(
        internal_model, to_coordinate_system, to_selector_id
    )

    crossmap = crossmap_to_hgvs_setup(to_coordinate_system, selector_model, degenerate)

    for point, path in yield_point_locations_for_main_reference(internal_model):
        set_by_path(hgvs_model, path, point_to_hgvs(point, **crossmap))

    return hgvs_model
