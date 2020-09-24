import copy

from mutalyzer_crossmapper import Coding, Genomic, NonCoding

from ..description import (
    get_reference_id,
    get_selector_id,
    yield_inserted_other_reference,
    yield_point_locations_for_main_reference,
)
from ..reference import get_coordinate_system_from_selector_id, get_selector_model
from ..util import set_by_path


def add_msg(dictionary, message_type, message):
    if dictionary.get(message_type) is None:
        dictionary[message_type] = []
    dictionary[message_type].append(message)


def check_selector_model(description_model, selector_model):
    if description_model["coordinate_system"] in ["c", "n"]:
        if not selector_model.get("exon"):
            message = {
                "code": "ESELECTORMODELNOEXONS",
                "detail": "No exons found in the reference for {} selector.".format(
                    selector_model["id"]
                ),
            }
            if description_model.get("reference"):
                add_msg(description_model["reference"]["selector"], "errors", message)
            elif description_model.get("source"):
                add_msg(description_model["source"]["selector"], "errors", message)

    if description_model["coordinate_system"] in ["c"]:
        if not selector_model.get("cds"):
            message = {
                "code": "ESELECTORMODELNOCDS",
                "detail": "No CDS found in the reference for {} selector.".format(
                    selector_model["id"]
                ),
            }
            if description_model.get("reference"):
                add_msg(description_model["reference"]["selector"], "errors", message)
            elif description_model.get("source"):
                add_msg(description_model["source"]["selector"], "errors", message)


def get_point_value(point):
    value = point["position"]
    if point.get("outside_cds") and point["outside_cds"] == "upstream":
        value *= -1
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


def point_to_coding(point, crossmap_function, point_function):
    if point.get("uncertain"):
        return {"type": "point", "uncertain": True}
    else:
        return create_exact_point_model(crossmap_function(point_function(point)))


def point_to_internal(point, crossmap):
    new_point = point_to_coding(point, **crossmap)
    return new_point


def crossmap_to_internal_setup(coordinate_system, selector_model=None):
    if coordinate_system == "g":
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
        }
    elif coordinate_system == "n":
        crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
        return {
            "crossmap_function": crossmap.noncoding_to_coordinate,
            "point_function": point_to_x_coding,
        }


def initialize_internal_model(description_model):
    internal_description_model = copy.deepcopy(description_model)
    internal_description_model["coordinate_system"] = "x"
    return internal_description_model


def points_to_internal_coordinates(description, references):
    reference_id = get_reference_id(description)
    coordinate_system = description.get("coordinate_system")
    selector_id = get_selector_id(description)
    selector_model = (
        get_selector_model(references[reference_id]["annotations"], selector_id)
        if selector_id
        else None
    )

    internal_description_model = initialize_internal_model(description)
    crossmap = crossmap_to_internal_setup(coordinate_system, selector_model)

    for point, path in yield_point_locations_for_main_reference(description):
        set_by_path(
            internal_description_model, path, point_to_internal(point, crossmap)
        )

    return internal_description_model


def to_internal_coordinates(description, references):
    """

    :param description: Description model.
    :param references: Dictionary with reference models with their ids as keys.
    :return: Converted description model with locations in the internal coordinate system.
    """
    internal_model = points_to_internal_coordinates(description, references)
    for inserted, path in yield_inserted_other_reference(description):
        set_by_path(
            internal_model, path, points_to_internal_coordinates(inserted, references)
        )

    return internal_model
