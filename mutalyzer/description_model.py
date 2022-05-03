def get_reference_id(model):
    """
    Get the reference ID from a description model (inserted can be
    considered also as a description model input).

    :param model: Description model (inserted works also).
    :return: The ID of the reference, if found, otherwise None.
    """
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
    :return: The ID of the selector, if found, otherwise None.
    """
    if (
        model.get("reference")
        and model["reference"].get("selector")
        and model["reference"]["selector"].get("id")
    ):
        return model["reference"]["selector"]["id"]
    elif (
        model.get("source")
        and isinstance(model["source"], dict)
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


def yield_ranges_main_reference(model, path=[]):
    for k in model.keys():
        if k in ["location", "start", "end"]:
            if model[k]["type"] == "range":
                yield model[k], path + [k]
            else:
                yield from yield_ranges_main_reference(model[k], path + [k])
        elif k in ["variants", "inserted"]:
            for i, sub_model in enumerate(model[k]):
                if not isinstance(sub_model.get("source"), dict):
                    yield from yield_ranges_main_reference(sub_model, path + [k, i])


def yield_point_locations_all(model, path=[]):
    for k in model.keys():
        if k in ["location", "start", "end"]:
            if model[k]["type"] == "point":
                yield model[k], path + [k]
            else:
                yield from yield_point_locations_all(model[k], path + [k])
        elif k in ["variants", "inserted", "deleted"]:
            for i, sub_model in enumerate(model[k]):
                yield from yield_point_locations_all(sub_model, path + [k, i])


def yield_sub_model(model, keys, types=None, path=[]):
    """

    :param model:
    :param keys:
    :param types:
    :param path:
    """
    if isinstance(model, dict):
        for k in model.keys():
            if (k in keys and not types) or (
                k in keys and types and model[k]["type"] in types
            ):
                yield model[k], path + [k]
            yield from yield_sub_model(model[k], keys, types, path + [k])
    elif isinstance(model, list):
        for i, sub_model in enumerate(model):
            yield from yield_sub_model(sub_model, keys, types, path + [i])


def yield_values(model, keys, path=[]):
    if isinstance(model, dict):
        for k in model.keys():
            if k in keys:
                yield model[k], path + [k]
            yield from yield_values(model[k], keys, path + [k])
    elif isinstance(model, list):
        for i, sub_model in enumerate(model):
            yield from yield_values(sub_model, keys, path + [i])


def get_locations_min_max(model):
    """
    Get the minimum and maximum positions from all the locations present in
    the model.
    :param model: Desription model.
    :return: Minimum and maximum positions.
    """
    locations = [
        x[0]["position"]
        for x in yield_point_locations_for_main_reference_variants(model)
        if x[0].get("position")
    ]
    if locations:
        return min(locations), max(locations)
    else:
        return None, None


def model_to_string(model, exclude_superfluous_selector=True, aa="verbatim"):
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
    else:
        reference_id = None
    if reference_id:
        selector_id = get_selector_id(model)
        if (
            reference_id
            and exclude_superfluous_selector
            and reference_id == selector_id
        ):
            selector_id = None
        if selector_id:
            reference = "{}({})".format(reference_id, selector_id)
        else:
            reference = "{}".format(reference_id)
    else:
        reference = None
    if model.get("coordinate_system"):
        coordinate_system = model.get("coordinate_system") + "."
    else:
        coordinate_system = ""
    if isinstance(model.get("variants"), list):
        if model.get("type") == "description_protein":
            variants = variants_to_description(model["variants"], True)
        else:
            variants = variants_to_description(model["variants"], False)
        if model.get("predicted"):
            variants = f"({variants})"
        if reference:
            return "{}:{}{}".format(reference, coordinate_system, variants)
        else:
            return "{}".format(variants)
    if model.get("location"):
        return "{}:{}{}".format(
            reference, coordinate_system, location_to_description(model.get("location"))
        )


def variants_to_description(variants, protein=False, aa="verbatim"):
    """
    Convert a list of variant models to string.
    :param variants: Variant models.
    :return: Variants string representation.
    """
    if isinstance(variants, list):
        if len(variants) == 0:
            return "="
        variants_list = []
        for variant in variants:
            variants_list.append(variant_to_description(variant, protein))
        if len(variants_list) > 1:
            return "[{}]".format(";".join(variants_list))
        elif len(variants_list) == 1:
            return variants_list[0]


def variant_to_description(variant, protein=False, aa="verbatim"):
    """
    Convert the variant dictionary model to string.
    :param variant: Variant model.
    :return: Variant model string representation.
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

    if variant.get("deleted"):
        deleted = inserted_to_description(variant["deleted"])

    variant_type = variant.get("type")
    if variant_type == "substitution":
        variant_type = deleted
        if not protein:
            variant_type += ">"
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


def inserted_to_description(inserted, aa="verbatim"):
    """
    Convert the inserted dictionary model to string.
    :param inserted: Inserted dictionary model.
    :return: Inserted string representation.
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


def location_to_description(location, aa="verbatim"):
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
            return "{}_{}".format(
                location_to_description(location.get("start")),
                location_to_description(location.get("end")),
            )


def point_to_description(point, aa="verbatim"):
    """
    Convert the position dictionary model to string.
    :param point: Position dictionary.
    :return: Equivalent position string representation.
    """
    outside_cds = offset = ""
    if point.get("amino_acid"):
        sequence = point.get("amino_acid")
    else:
        sequence = ""
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
    return "{}{}{}{}".format(sequence, outside_cds, position, offset)


def length_to_description(length):
    """
    Convert the length dictionary model to string.
    :param length: Length dictionary model.
    :return: Equivalent length string representation.
    """
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
