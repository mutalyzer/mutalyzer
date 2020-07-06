import copy
from mutalyzer_crossmapper import Crossmap

from normalizer.util import get_start, get_end
from normalizer.reference import get_mol_type, get_selector_model


def range_to_point(range_location):
    """
    Convert a point location to a range location.

    :param range_location: A point location model complying object.
    :return: A point location model complying object.
    """
    if get_start(range_location) == get_end(range_location):
        return {"type": "point", "position": get_start(range_location)}
    else:
        return range_location


def location_to_hgvs(location, crossmap_function, point_function):
    new_location = copy.deepcopy(location)
    if location["type"] == "range":
        new_location["start"] = point_function(
            crossmap_function(location["start"]["position"])
        )
        new_location["end"] = point_function(
            crossmap_function(location["end"]["position"])
        )
    elif location["type"] == "point":
        new_location = point_function(crossmap_function(location["position"]))
    return new_location


def location_to_hgvs_ins_fix(location, variant_type):
    new_location = copy.deepcopy(location)

    if variant_type == "insertion":
        new_location["start"]["position"] -= 1
    elif location["type"] == "range":
        new_location["end"]["position"] -= 1
        if get_start(new_location) == get_end(new_location):
            new_location = range_to_point(new_location)

    return new_location


def genomic_to_point(genomic):
    point = {"type": "point", "position": genomic}
    return point


def crossmap_coordinate_to_genomic_setup():
    crossmap = Crossmap()
    return {
        "crossmap_function": crossmap.coordinate_to_genomic,
        "point_function": genomic_to_point,
    }


def coordinate_to_genomic(variants):
    crossmap = crossmap_coordinate_to_genomic_setup()
    new_variants = []

    for variant in variants:
        new_variant = copy.deepcopy(variant)
        hgvs_indexing_location = location_to_hgvs_ins_fix(
            variant["location"], variant["type"]
        )
        new_variant["location"] = location_to_hgvs(hgvs_indexing_location, **crossmap)
        new_variants.append(new_variant)

    return new_variants


def fix_crossmap_coding(exon, coding, crossmap):
    first_exon_coding = crossmap.coordinate_to_coding(exon)
    if coding[0] == first_exon_coding[0]:
        coding = (coding[0] + coding[1], 0, coding[2])
    return coding


def coordinate_to_coding(variants, selector_model):
    crossmap = Crossmap(
        selector_model["exon"],
        selector_model["cds"][0],
        inverted=selector_model["inverted"],
    )
    coordinate_variants = []
    for variant in variants:
        new_variant = copy.deepcopy(variant)
        variant["location"] = location_to_hgvs_ins_fix(
            variant["location"], variant["type"]
        )
        if variant["location"]["type"] == "point":
            coding = crossmap.coordinate_to_coding(variant["location"]["position"])
            # compensate for the crossmaper
            coding = fix_crossmap_coding(selector_model["exon"][0][0], coding, crossmap)
            coding = fix_crossmap_coding(
                selector_model["exon"][-1][1], coding, crossmap
            )
            new_variant["location"] = coding_to_point(coding)
        if variant["location"]["type"] == "range":
            # compensate for the crossmaper
            first_exon_coding = crossmap.coordinate_to_coding(
                selector_model["exon"][0][0]
            )
            coding_start = crossmap.coordinate_to_coding(
                variant["location"]["start"]["position"]
            )
            if coding_start[0] == first_exon_coding[0]:
                coding_start = (coding_start[0] + coding_start[1], 0, 0)
            coding_end = crossmap.coordinate_to_coding(
                variant["location"]["end"]["position"]
            )
            if coding_end[0] == first_exon_coding[0]:
                coding_end = (coding_end[0] + coding_end[1], 0, 0)

            new_variant["location"] = {
                "type": "range",
                "start": coding_to_point(coding_start),
                "end": coding_to_point(coding_end),
            }
        coordinate_variants.append(new_variant)
    return coordinate_variants


def identify_coordinate_system(selector_model):
    if selector_model["type"] in ["mRNA"]:
        return "c"
    if selector_model["type"] in ["ncRNA"]:
        return "n"


def coding_to_point(coding):
    position, offset, section = coding
    point = {"type": "point", "position": position}

    if section == 0:
        point["outside_cds"] = "upstream"
        point["position"] *= -1
    elif section == 2:
        point["outside_cds"] = "downstream"

    if offset != 0:
        point["offset"] = {"value": offset}

    return point


def to_hgvs_locations(variants, reference_model, selector_id=None):
    """
    Converts the variant locations present from the internal coordinate system
    to the hgvs locations, according the the selector_id type.

    :param variants: Variants with locations in the crossmap coordinate system.
    :param reference_model: Reference models, dictionary.
    :param selector_id: The selector id.
    :return: Variants with locations in selector's HGVS coordinate system.
    """
    if selector_id is None and get_mol_type(reference_model):
        coordinate_system = "g"
    else:
        selector_model = get_selector_model(reference_model["model"], selector_id)
        coordinate_system = identify_coordinate_system(selector_model)

    if coordinate_system is None:
        raise Exception("No coordinate system. Unable to convert.")

    if coordinate_system == "g":
        hgvs_variants = coordinate_to_genomic(variants)
    elif coordinate_system == "c":
        hgvs_variants = coordinate_to_coding(variants, selector_model)
    elif coordinate_system == "n":
        selector_model["cds"] = [
            (selector_model["exon"][0][0], selector_model["exon"][-1][-1])
        ]
        hgvs_variants = coordinate_to_coding(variants, selector_model)
    else:
        raise Exception("Unsupported coordinate system: {}.".format(coordinate_system))

    if selector_id is None:
        reference = {"id": reference_model["model"]["id"]}
    else:
        reference = {
            "id": reference_model["model"]["id"],
            "selector": {"id": selector_id},
        }
    return {
        "reference": reference,
        "variants": hgvs_variants,
        "coordinate_system": coordinate_system,
    }
