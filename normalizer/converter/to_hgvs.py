import copy
from mutalyzer_crossmapper import Genomic, NonCoding, Coding

from normalizer.util import get_start, get_end, create_exact_point_model
from normalizer.reference import get_mol_type, get_selector_model


def range_to_point(range_location):
    """
    Convert a point location to a range location.

    :param range_location: A point location model complying object.
    :return: A point location model complying object.
    """
    if range_location.get("uncertain") or range_location["start"].get("uncertain") or range_location["end"].get("uncertain"):
        return range_location
    elif range_location["start"] == range_location["end"]:
        return range_location["start"]
    return range_location


def genomic_to_point(genomic):
    point = {"type": "point", "position": genomic}
    return point


def identify_coordinate_system(selector_model):
    if selector_model["type"] in ["mRNA"]:
        return "c"
    if selector_model["type"] in ["ncRNA"]:
        return "n"


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
    print(coding)
    print(point)
    return point


def noncoding_to_point(noncoding):
    position, offset, section = noncoding[:3]
    point = {"type": "point", "position": position}

    if offset != 0:
        point["offset"] = {"value": offset}
    return point


def crossmap_coordinate_to_genomic_setup():
    crossmap = Genomic()
    return {
        "crossmap_function": crossmap.coordinate_to_genomic,
        "point_function": genomic_to_point,
    }


def crossmap_coordinate_to_coding_setup(selector, degenerate=False):
    crossmap = Coding(selector["exon"], selector["cds"][0], selector["inverted"])
    return {
        "crossmap_function": crossmap.coordinate_to_coding,
        "point_function": coding_to_point,
        "degenerate": degenerate
    }


def crossmap_coordinate_to_noncoding_setup(selector, degenerate=False):
    crossmap = NonCoding(selector["exon"], selector["inverted"])
    return {
        "crossmap_function": crossmap.coordinate_to_noncoding,
        "point_function": noncoding_to_point,
        "degenerate": degenerate
    }


def crossmap_to_hgvs_setup(coordinate_system, selector_model=None, degenerate=False):
    """
    Returns a crossmap instance able to convert from the internal system
    to the to hgvs system.
    """
    if coordinate_system == "g":
        crossmap = crossmap_coordinate_to_genomic_setup()
    elif coordinate_system == "c":
        crossmap = crossmap_coordinate_to_coding_setup(selector_model, degenerate)
    elif coordinate_system == "n":
        crossmap = crossmap_coordinate_to_noncoding_setup(selector_model)
    else:
        raise Exception("Unsupported coordinate system: {}.".format(coordinate_system))

    return crossmap


def inserted_to_hgvs(inserted):
    crossmap = {
        "crossmap_function": Genomic().coordinate_to_genomic,
        "point_function": genomic_to_point,
    }
    return location_to_hgvs(
        location=inserted["location"], variant_type=None, crossmap=crossmap
    )


def update_interval(range_location, variant_type):
    if (
        variant_type == "insertion"
        and range_location["start"]["type"] == "point"
        and not range_location["start"].get("uncertain")
    ):
        range_location["start"]["position"] -= 1
    elif (
        not range_location["end"].get("uncertain")
        and range_location["end"]["type"] == "point"
    ):
        range_location["end"]["position"] -= 1


def point_to_hgvs(point, crossmap_function, point_function, degenerate=False):
    if point.get("uncertain"):
        return {"type": "point", "uncertain": True}
    else:
        if degenerate:
            return point_function(crossmap_function(point["position"], degenerate))
        else:
            return point_function(crossmap_function(point["position"]))


def range_to_hgvs(range_location, variant_type, crossmap):
    update_interval(range_location, variant_type)
    new_range = {
        "start": point_to_hgvs(range_location["start"], **crossmap),
        "end": point_to_hgvs(range_location["end"], **crossmap),
        "type": "range",
    }
    if range_location.get("uncertain"):
        new_range["uncertain"] = True
    return new_range


def location_to_hgvs(location, variant_type, crossmap):
    new_location = copy.deepcopy(location)
    if location["start"]["type"] == "point" and location["end"]["type"] == "point":
        new_location = range_to_hgvs(new_location, variant_type, crossmap)
    else:
        update_interval(new_location, variant_type)
        if location["start"]["type"] == "range":
            new_location["start"] = range_to_hgvs(new_location["start"], variant_type, crossmap)
        else:
            new_location["start"] = point_to_hgvs(new_location["start"], **crossmap)
        if location["end"]["type"] == "range":
            new_location["end"] = range_to_hgvs(new_location["end"], variant_type, crossmap)
        else:
            new_location["end"] = point_to_hgvs(new_location["end"], **crossmap)
    new_location = range_to_point(new_location)
    return new_location


def variants_to_hgvs(variants, crossmap):
    new_variants = []

    for variant in variants:
        new_variant = copy.deepcopy(variant)
        new_variant["location"] = location_to_hgvs(
            variant["location"], variant['type'], crossmap
        )
        if new_variant.get("inserted"):
            for ins in new_variant["inserted"]:
                if ins.get("location"):
                    ins["location"] = inserted_to_hgvs(ins)
        new_variants.append(new_variant)

    return new_variants


def to_hgvs_locations(variants, reference_model, selector_id=None, degenerate=False):
    """
    Converts the variant locations present from the internal coordinate system
    to the hgvs locations, according the the selector_id type.

    :param variants: Variants with locations in the crossmap coordinate system.
    :param reference_model: Reference models, dictionary.
    :param selector_id: The selector id.
    :param degenerate: Whether to output degenerate positions in c. or n.
    :return: Variants with locations in selector's HGVS coordinate system.
    """
    selector_model = None
    if selector_id is None and get_mol_type(reference_model):
        coordinate_system = "g"
    else:
        selector_model = get_selector_model(reference_model["model"], selector_id)
        coordinate_system = identify_coordinate_system(selector_model)

    if coordinate_system is None:
        raise Exception("No coordinate system. Unable to convert.")

    crossmap = crossmap_to_hgvs_setup(coordinate_system, selector_model, degenerate)

    if selector_id is None:
        reference = {"id": reference_model["model"]["id"]}
    else:
        reference = {
            "id": reference_model["model"]["id"],
            "selector": {"id": selector_id},
        }
    return {
        "reference": reference,
        "variants": variants_to_hgvs(variants, crossmap),
        "coordinate_system": coordinate_system,
    }
