import copy

from .util import get_end, get_start, set_by_path


def get_reference_id(model):
    if model.get("reference") and model["reference"].get("id"):
        return model["reference"]["id"]
    elif (
        model.get("source")
        and isinstance(model["source"], dict)
        and model["source"].get("id")
    ):
        return model["source"]["id"]


def get_selector_id(model):
    """
    Get the main selector ID from the description model (inserted can be
    considered also as a description model input). At the moment, no nesting,
    i.e., selector(selector(...)), is supported.

    :param model: Description model (inserted works also).
    :return: The ID of the selector, if provided, otherwise None.
    """
    if (
        model.get("reference")
        and model["reference"].get("selector")
        and model["reference"]["selector"].get("id")
    ):
        return model["reference"]["selector"]["id"]
    elif (
        model.get("source")
        and model["source"].get("selector")
        and model["source"]["selector"].get("id")
    ):
        return model["source"]["selector"]["id"]


def yield_reference_ids(model, path=[]):
    for k in model.keys():
        if k in ["reference", "source"]:
            if isinstance(model[k], dict) and model[k].get("id"):
                yield model[k]["id"], tuple(path + [k, "id"])
        elif k in ["variants", "inserted"]:
            for i, sub_model in enumerate(model[k]):
                yield from yield_reference_ids(sub_model, path + [k, i])


def yield_inserted_other_reference(model, path=[]):
    for k in model.keys():
        if k == "variants":
            for i, variant in enumerate(model[k]):
                yield from yield_inserted_other_reference(variant, path + [k, i])
        elif k == "inserted":
            for i, inserted in enumerate(model[k]):
                if isinstance(inserted.get("source"), dict):
                    yield inserted, tuple(path + [k, i])


def yield_reference_selector_ids(model, path=[]):
    for k in model.keys():
        if k in ["reference", "source"]:
            if (
                isinstance(model[k], dict)
                and model[k].get("id")
                and model[k].get("selector")
                and model[k]["selector"].get("id")
            ):
                yield (
                    model[k]["id"],
                    model[k]["selector"]["id"],
                    tuple(path + [k, "selector", "id"]),
                )
        elif k in ["variants", "inserted"]:
            for i, sub_model in enumerate(model[k]):
                yield from yield_reference_selector_ids(sub_model, path + [k, i])


def yield_reference_selector_ids_coordinate_system(model, path=[]):
    for k in model.keys():
        if k in ["reference", "source"] and isinstance(model[k], dict):
            c_s = model.get("coordinate_system")
            c_s_path = tuple(path + ["coordinate_system"])
            r_id = model[k].get("id")
            r_path = tuple(path + [k, "id"])
            s_id = None
            s_path = tuple(path + [k])
            if model[k].get("selector") and model[k]["selector"].get("id"):
                s_id = model[k]["selector"]["id"]
                s_path = tuple(path + [k, "selector", "id"])
            yield c_s, c_s_path, r_id, r_path, s_id, s_path
        elif k in ["variants", "inserted"]:
            for i, sub_model in enumerate(model[k]):
                yield from yield_reference_selector_ids_coordinate_system(
                    sub_model, path + [k, i]
                )


def yield_point_locations_for_main_reference(model, path=[]):
    for k in model.keys():
        if k in ["location", "start", "end"]:
            if model[k]["type"] == "point":
                yield model[k], path + [k]
            else:
                yield from yield_point_locations_for_main_reference(
                    model[k], path + [k]
                )
        elif k in ["variants", "inserted"]:
            for i, sub_model in enumerate(model[k]):
                if not isinstance(sub_model.get("source"), dict):
                    yield from yield_point_locations_for_main_reference(
                        sub_model, path + [k, i]
                    )


def yield_point_locations_for_main_reference_variants(model, path=[]):
    for k in model.keys():
        if k in ["location", "start", "end"]:
            if model[k]["type"] == "point":
                yield model[k], path + [k]
            else:
                yield from yield_point_locations_for_main_reference(
                    model[k], path + [k]
                )
        elif k in ["variants"]:
            for i, sub_model in enumerate(model[k]):
                if not isinstance(sub_model.get("source"), dict):
                    yield from yield_point_locations_for_main_reference_variants(
                        sub_model, path + [k, i]
                    )


def yield_range_locations_for_main_reference(model, path=[]):
    for k in model.keys():
        if k in ["location", "start", "end"]:
            if model[k]["type"] == "range":
                yield model[k], path + [k]
            else:
                yield from yield_range_locations_for_main_reference(
                    model[k], path + [k]
                )
        elif k in ["variants", "inserted"]:
            for i, sub_model in enumerate(model[k]):
                if not isinstance(sub_model.get("source"), dict):
                    yield from yield_range_locations_for_main_reference(
                        sub_model, path + [k, i]
                    )


def yield_view_nodes(model, path=[]):
    for k in model.keys():
        if k in ["reference", "location", "start", "end", "selector"]:
            yield from yield_view_nodes(model[k], path + [k])
        elif k in ["variants", "inserted", "deleted"]:
            for i, sub_model in enumerate(model[k]):
                yield from yield_view_nodes(sub_model, path + [k, i])
        elif k == "source" and isinstance(model[k], dict):
            yield from yield_view_nodes(model[k], path + [k])
        if k in ["position", "id", "coordinate_system"]:
            yield model[k], path + [k]


def get_view_model(model):
    view_model = copy.deepcopy(model)
    for view_value, path in yield_view_nodes(model):
        set_by_path(view_model, path, {"view": view_value})

    return view_model


def model_to_string(model, exclude_superfluous_selector=True):
    """
    Convert the variant description model to string.
    :param exclude_superfluous_selector: Do not include the selector_id
           if is it the same as the reference_id.
    :param model: Dictionary holding the variant description model.
    :return: Equivalent reference string representation.
    """
    if model.get("reference"):
        reference_id = model["reference"]["id"]
    elif model.get("source"):
        reference_id = model["source"]["id"]
    selector_id = get_selector_id(model)
    if exclude_superfluous_selector and reference_id == selector_id:
        selector_id = None
    if selector_id:
        reference = "{}({})".format(reference_id, selector_id)
    else:
        reference = "{}".format(reference_id)
    if model.get("coordinate_system"):
        coordinate_system = model.get("coordinate_system") + "."
    else:
        coordinate_system = ""
    if model.get("variants"):
        return "{}:{}{}".format(
            reference, coordinate_system, variants_to_description(model.get("variants"))
        )
    if model.get("location"):
        return "{}:{}{}".format(
            reference, coordinate_system, location_to_description(model.get("location"))
        )


def reference_to_description(reference):
    """
    Convert the reference dictionary model to string.
    :param reference: Dictionary holding the reference model.
    :return: Equivalent reference string representation.
    """
    version = ""
    if isinstance(reference, dict):
        if reference.get("type") == "genbank":
            accession = reference.get("accession")
            if reference.get("version"):
                version = ".{}".format(reference["version"])
        elif reference.get("type") == "lrg":
            accession = reference.get("id")
    return "{}{}".format(accession, version)


def variants_to_description(variants, sequences=None):
    if isinstance(variants, list):
        variants_list = []
        for variant in variants:
            variants_list.append(variant_to_description(variant, sequences))
        if len(variants_list) > 1:
            return "[{}]".format(";".join(variants_list))
        elif len(variants_list) == 1:
            return variants_list[0]


def variant_to_description(variant, sequences=None):
    """
    Convert the variant dictionary model to string.
    :return: Equivalent variant string representation.
    """
    deleted_location = deleted = inserted = ""
    if variant.get("location"):
        deleted_location = location_to_description(variant.get("location"))
    if variant.get("inserted"):
        inserted = inserted_to_description(variant["inserted"])
        if (
            variant["type"] == "repeat"
            and len(variant.get("inserted")) == 1
            and (
                (
                    variant["inserted"][0].get("location")
                    and not variant["inserted"][0].get("sequence")
                )
                or (
                    variant["inserted"][0].get("length")
                    and variant["inserted"][0]["length"].get("value")
                )
            )
        ):
            inserted = "[{}]".format(inserted)
    variant_type = variant.get("type")

    if variant.get("deleted"):
        deleted = inserted_to_description(variant["deleted"])

    if variant_type == "substitution":
        variant_type = deleted + ">"
    elif variant_type == "deletion":
        variant_type = "del" + deleted
    elif variant_type == "deletion_insertion":
        variant_type = "del" + deleted + "ins"
    elif variant_type == "insertion":
        variant_type = "ins"
    elif variant_type == "duplication":
        variant_type = "dup" + inserted
        inserted = ""
    elif variant_type == "inversion":
        variant_type = "inv"
    elif variant_type == "conversion":
        variant_type = "con"
    elif variant_type == "equal":
        variant_type = "="
    else:
        variant_type = ""
    return "{}{}{}".format(deleted_location, variant_type, inserted)


def inserted_to_description(inserted):
    """
    Convert the insertions dictionary model to string.
    :param inserted: Insertions dictionary.
    :return: Equivalent insertions string representation.
    """
    descriptions = []
    for insert in inserted:
        if insert.get("sequence"):
            descriptions.append(insert["sequence"])
        elif insert.get("source") and isinstance(insert["source"], dict):
            descriptions.append(model_to_string(insert))
        elif insert.get("location"):
            descriptions.append(location_to_description(insert["location"]))
        elif insert.get("length"):
            descriptions.append(length_to_description(insert["length"]))
        if insert.get("repeat_number"):
            descriptions[-1] += "[{}]".format(insert["repeat_number"]["value"])
        if insert.get("inverted"):
            descriptions[-1] += "inv"
    if len(inserted) > 1:
        return "[{}]".format(";".join(descriptions))
    else:
        return descriptions[0]


def location_to_description(location):
    """
    Convert the location dictionary model to string.
    :param location: Location dictionary.
    :return: Equivalent location string representation.
    """
    if location["type"] == "point":
        return point_to_description(location)
    if location["type"] == "range":
        if location.get("uncertain"):
            return "({}_{})".format(
                point_to_description(location.get("start")),
                point_to_description(location.get("end")),
            )
        else:
            start = location_to_description(location.get("start"))
            end = location_to_description(location.get("end"))
            return "{}_{}".format(start, end)


def point_to_description(point):
    """
    Convert the position dictionary model to string.
    :param point: Position dictionary.
    :return: Equivalent position string representation.
    """
    outside_cds = offset = ""
    if point.get("outside_cds"):
        if point["outside_cds"] == "downstream":
            outside_cds = "*"
        elif point["outside_cds"] == "upstream":
            outside_cds = "-"
    if point.get("uncertain"):
        position = "?"
    else:
        position = str(point.get("position"))
    if point.get("offset"):
        if point["offset"].get("value"):
            offset = "%+d" % point["offset"]["value"]
        elif point["offset"].get("uncertain"):
            if point["offset"].get("upstream"):
                offset = "-?"
            elif point["offset"].get("downstream"):
                offset = "+?"
    return "{}{}{}".format(outside_cds, position, offset)


def length_to_description(length):
    if length["type"] == "point":
        if length.get("value"):
            return str(length["value"])
        elif length.get("uncertain"):
            return "?"
    if length["type"] == "range":
        output = "{}_{}".format(
            length_to_description(length.get("start")),
            length_to_description(length.get("end")),
        )
        if length.get("uncertain"):
            return "({})".format(output)
        else:
            return output


def get_locations_start_end(model):
    locations = [
        x[0]["position"]
        for x in yield_point_locations_for_main_reference_variants(model)
        if x[0].get("position")
    ]
    return min(locations), max(locations)
