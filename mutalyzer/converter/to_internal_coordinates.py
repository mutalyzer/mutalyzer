import copy

from mutalyzer_crossmapper import Coding, Genomic, NonCoding

from ..description_model import (
    get_reference_id,
    get_selector_id,
    yield_inserted_other_reference,
    yield_point_locations_for_main_reference,
    yield_ranges_main_reference,
)
from ..reference import (
    get_coordinate_system_from_reference,
    get_coordinate_system_from_selector_id,
    get_internal_selector_model,
)
from ..util import set_by_path


def get_point_value(point):
    value = point["position"]
    if point.get("outside_cds") and point["outside_cds"] == "upstream":
        value *= -1
    if point.get("offset") and point["offset"].get("value"):
        value += point["offset"]["value"]
    return value


def point_to_x_coding(point):
    position = point["position"]
    if point.get("outside_cds"):
        if point["outside_cds"] == "upstream":
            section = -1
            position = -1 * position
        elif point["outside_cds"] == "downstream":
            section = 1
    else:
        section = 0
    if point.get("offset"):
        offset = point["offset"]["value"]
    else:
        offset = 0
    return position, offset, section, 0


def create_exact_point_model(point):
    return {"type": "point", "position": point}


def point_to_coding(point, crossmap_function, point_function, inverted=False):
    if point.get("uncertain"):
        return {"type": "point", "uncertain": True}
    else:
        internal_point = crossmap_function(point_function(point))
        if inverted and point.get("shift"):
            internal_point += point["shift"]
        point_model = {"shift": point["shift"]} if point.get("shift") else {}
        point_model.update(create_exact_point_model(internal_point))
        return point_model


def point_to_internal(point, crossmap):
    new_point = point_to_coding(point, **crossmap)
    if point.get("amino_acid"):
        new_point["amino_acid"] = point["amino_acid"]
    return new_point


def crossmap_to_internal_setup(coordinate_system, selector_model=None):
    if coordinate_system in ["g", "m", "p", None]:
        crossmap = Genomic()
        return {
            "crossmap_function": crossmap.genomic_to_coordinate,
            "point_function": get_point_value,
        }
    elif coordinate_system == "c":
        crossmap = Coding(
            selector_model["exon"],
            selector_model["cds"][0],
            selector_model["inverted"],
        )
        return {
            "crossmap_function": crossmap.coding_to_coordinate,
            "point_function": point_to_x_coding,
            "inverted": selector_model["inverted"],
        }
    elif coordinate_system == "n":
        crossmap = Coding(
            selector_model["exon"],
            (selector_model["exon"][0][0], selector_model["exon"][-1][-1]),
            selector_model["inverted"],
        )
        return {
            "crossmap_function": crossmap.coding_to_coordinate,
            "point_function": point_to_x_coding,
            "inverted": selector_model["inverted"],
        }


def initialize_internal_model(model):
    internal_model = copy.deepcopy(model)
    if internal_model.get("reference"):
        internal_model["reference"] = {"id": get_reference_id(model)}
    if internal_model.get("source") and internal_model["source"].get("selector"):
        internal_model["source"].pop("selector")
    internal_model["coordinate_system"] = "x"
    return internal_model


def invert_sequences(variant, element_type):
    if variant.get(element_type):
        variant[element_type].reverse()
        for element in variant[element_type]:
            if element.get("sequence"):
                if element.get("inverted"):
                    element.pop("inverted")
                else:
                    element["inverted"] = True


def reverse_strand(internal_model):
    internal_model["variants"].reverse()
    for variant in internal_model["variants"]:
        invert_sequences(variant, "deleted")
        invert_sequences(variant, "inserted")

    for loc, path in yield_ranges_main_reference(internal_model):
        loc["start"], loc["end"] = loc["end"], loc["start"]


def get_coordinate_system(model, references):
    coordinate_system = model.get("coordinate_system")
    if coordinate_system == "r":
        if get_selector_id(model):
            coordinate_system = get_coordinate_system_from_selector_id(
                references[get_reference_id(model)], get_selector_id(model)
            )
        else:
            coordinate_system = get_coordinate_system_from_reference(
                references[get_reference_id(model)]
            )
    return coordinate_system


def points_to_internal_coordinates(model, references):
    reference_id = get_reference_id(model)
    coordinate_system = get_coordinate_system(model, references)
    selector_id = get_selector_id(model)
    selector_model = (
        get_internal_selector_model(
            references[reference_id]["annotations"], selector_id, True
        )
        if selector_id
        else None
    )

    internal_model = initialize_internal_model(model)
    crossmap = crossmap_to_internal_setup(coordinate_system, selector_model)

    for point, path in yield_point_locations_for_main_reference(model):
        set_by_path(internal_model, path, point_to_internal(point, crossmap))

    return internal_model


def to_internal_coordinates(model, references):
    """

    :param model: Description model.
    :param references: Dictionary with reference models with their ids as keys.
    :return: Converted description model with locations in the internal coordinate system.
    """
    internal_model = points_to_internal_coordinates(model, references)
    for inserted, path in yield_inserted_other_reference(model):
        set_by_path(
            internal_model, path, points_to_internal_coordinates(inserted, references)
        )

    reference_id = get_reference_id(model)
    selector_id = get_selector_id(model)
    if selector_id:
        selector_model = (
            get_internal_selector_model(
                references[reference_id]["annotations"], selector_id, True
            )
            if selector_id
            else None
        )

        if selector_model and selector_model.get("inverted"):
            reverse_strand(internal_model)

    return internal_model
