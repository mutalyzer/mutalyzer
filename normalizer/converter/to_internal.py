import copy
import warnings

from mutalyzer_crossmapper import Genomic, NonCoding, Coding

from normalizer.description import *
from normalizer.reference import get_available_selectors, get_selector_model
from normalizer.util import create_exact_point_model, sort_location_tuples


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


def fix_selector_id(reference_models, reference_id, coordinate_system):
    available_selectors = get_available_selectors(
        reference_models[reference_id]["model"], coordinate_system
    )
    if len(available_selectors) == 0:
        raise Exception(
            "ENOSELECTOR: {} coordinate system used but no selector ID "
            "provided in the description. In addition, there is no "
            "selector available in the reference model.".format(coordinate_system)
        )
    elif len(available_selectors) == 1:
        warnings.warn(
            "WNOSELECTOR: {} coordinate system used but no selector ID "
            "provided in the description. Only {} present in the reference,"
            " which is chosen as default.".format(
                coordinate_system, available_selectors[0]
            )
        )
        return available_selectors[0]
    elif len(available_selectors) > 1:
        raise Exception(
            "ENOSELECTOR: {} coordinate system used but no selector ID "
            "provided in the description. Please choose between the "
            "following selectors available in the reference: {}".format(
                coordinate_system, available_selectors
            )
        )


def get_point_value(point):
    value = point["position"]
    if point.get("outside_cds") and point["outside_cds"] == "upstream":
        value *= -1
    return value


def crossmap_genomic_to_coordinate_setup():
    crossmap = Genomic()
    return {
        "crossmap_function": crossmap.genomic_to_coordinate,
        "point_function": get_point_value,
    }


def crossmap_coding_to_coordinate_setup(description, references):
    reference_id = description["reference"]["id"]

    selector_id = get_selector_id(description)
    if selector_id is None:
        selector_id = fix_selector_id(references, reference_id, "c")
    selector = get_selector_model(references[reference_id]["model"], selector_id)
    # if not selector.get("exon") and

    crossmap = Coding(selector["exon"], selector["cds"][0], selector["inverted"])

    return {
        "crossmap_function": crossmap.coding_to_coordinate,
        "point_function": point_to_x_coding,
    }


def crossmap_noncoding_to_coordinate_setup(description, references):
    reference_id = description["reference"]["id"]

    selector_id = get_selector_id(description)
    if selector_id is None:
        selector_id = fix_selector_id(references, reference_id, "n")

    selector = get_selector_model(references[reference_id]["model"], selector_id)
    selector["exon"] = sort_location_tuples(selector["exon"])
    crossmap = NonCoding(selector["exon"], selector["inverted"])
    return {
        "crossmap_function": crossmap.noncoding_to_coordinate,
        "point_function": point_to_x_coding,
    }


def crossmap_to_x_setup(description, references):
    """
    Returns a crossmap instance able to convert from the coordinate system
    provided in the description model to the to internal system (crossmap
    coordinate).
    :param description: Description model.
    :param references: References models.
    :return: {'crossmap_function' : ..., 'point_function': ...}
    """
    coordinate_system = get_coordinate_system(description)
    if coordinate_system is None:
        warnings.warn("No coordinate system, we assume g.")
        # TODO: Improve message
        coordinate_system = "g"
        # TODO: Update also the description model.

    if coordinate_system == "g":
        crossmap = crossmap_genomic_to_coordinate_setup()
    elif coordinate_system == "c":
        crossmap = crossmap_coding_to_coordinate_setup(description, references)
    elif coordinate_system == "n":
        crossmap = crossmap_noncoding_to_coordinate_setup(description, references)
    else:
        raise Exception("Unsupported coordinate system: {}.".format(coordinate_system))

    return crossmap


def point_to_range(point_location):
    """
    Convert a point location to a range location.

    :param point_location: A point location model complying object.
    :return: A range location model complying object.
    """
    if point_location.get("uncertain"):
        return {
            "type": "range",
            "start": {"type": "point", "uncertain": True},
            "end": {"type": "point", "uncertain": True},
        }
    return {
        "type": "range",
        "start": {"type": "point", "position": point_location["position"]},
        "end": {"type": "point", "position": point_location["position"] + 1},
    }


def update_interval(range_location, variant_type):
    if (
        variant_type == "insertion"
        and range_location["start"]["type"] == "point"
        and not range_location["start"].get("uncertain")
    ):
        range_location["start"]["position"] += 1
    elif (
        not range_location["end"].get("uncertain")
        and range_location["end"]["type"] == "point"
    ):
        range_location["end"]["position"] += 1


def point_to_coding(point, crossmap_function, point_function):
    if point.get("uncertain"):
        return {"type": "point", "uncertain": True}
    else:
        return create_exact_point_model(crossmap_function(point_function(point)))


def point_to_internal(point, variant_type, crossmap):
    new_point = point_to_range(point_to_coding(point, **crossmap))
    return new_point


def range_to_internal(range_location, variant_type, crossmap):
    new_range = {
        "start": point_to_coding(range_location["start"], **crossmap),
        "end": point_to_coding(range_location["end"], **crossmap),
        "type": "range",
    }

    if range_location.get("uncertain"):
        new_range["uncertain"] = range_location["uncertain"]
    update_interval(new_range, variant_type)

    return new_range


def location_to_internal(location, variant_type, crossmap):
    """

    :param location:
    :param variant_type:
    :param crossmap:
    :return:
    """
    # import json
    # print(json.dumps(location, indent=2))
    if location["type"] == "range":
        if location["start"]["type"] == "range":
            new_location = {
                "start": range_to_internal(location["start"], variant_type, crossmap)
            }
        else:
            new_location = {"start": point_to_coding(location["start"], **crossmap)}
        if location["end"]["type"] == "range":
            new_location["end"] = range_to_internal(
                location["end"], variant_type, crossmap
            )
        else:
            new_location["end"] = point_to_coding(location["end"], **crossmap)
        new_location["type"] = "range"
        if location.get("uncertain"):
            new_location["uncertain"] = "uncertain"
        update_interval(new_location, variant_type)
    elif location["type"] == "point":
        new_location = point_to_internal(location, variant_type, crossmap)
    else:
        # Should never happen. TODO: Maybe raise an error?
        pass
    # print(json.dumps(new_location, indent=2))
    return new_location


def inserted_to_internal(inserted):
    crossmap = {
        "crossmap_function": Genomic().genomic_to_coordinate,
        "point_function": get_point_value,
    }
    return location_to_internal(
        location=inserted["location"], variant_type=None, crossmap=crossmap
    )


def variants_to_internal_locations(variants, crossmap):
    new_variants = []

    for variant in variants:
        new_variant = copy.deepcopy(variant)
        new_variant["location"] = location_to_internal(
            variant["location"], variant["type"], crossmap
        )
        if new_variant.get("inserted"):
            for ins in new_variant["inserted"]:
                if ins.get("location"):
                    ins["location"] = inserted_to_internal(ins)
        new_variants.append(new_variant)

    return new_variants


def to_internal_locations(description_model, references):
    """
    Converts the variant locations present in the description model to the
    internal coordinate system.

    :param description_model: Description model, dictionary.
    :param references: References models, dictionary.
    :return: Variants with locations in the internal coordinate system.
    """

    crossmap = crossmap_to_x_setup(description_model, references)

    return {
        "reference": {"id": description_model["reference"]["id"]},
        "coordinate_system": "x",
        "variants": variants_to_internal_locations(
            description_model["variants"], crossmap
        ),
    }
