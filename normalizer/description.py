import copy

from mutalyzer_hgvs_parser import parse_description_to_model

from .reference import Reference
from .util import add_msg


def get_reference_from_model(description_model):
    reference_id = get_reference_id(description_model)
    if not reference_id:
        return

    reference = Reference(reference_id)
    if not reference.model:
        if description_model.get("reference"):
            d_r = description_model["reference"]
        elif description_model.get("source"):
            d_r = description_model["source"]
        add_msg(
            d_r,
            "errors",
            {
                "code": "ERETR",
                "details": "Reference {} could not be retrieved.".format(reference_id),
            },
        )
    else:
        if reference.get_id() != reference_id:
            set_reference_id(description_model, reference.get_id())
        return reference


def get_references_from_description_model(model, references):
    if isinstance(model, dict):
        reference = get_reference_from_model(model)
        if reference:
            references[reference.get_id()] = reference
        for k in model.keys():
            if k in ["variants", "inserted"]:
                get_references_from_description_model(model[k], references)
    elif isinstance(model, list):
        for sub_model in model:
            get_references_from_description_model(sub_model, references)


def get_reference_id(model):
    if model.get("reference") and model["reference"].get("id"):
        return model["reference"]["id"]
    elif (
        model.get("source")
        and isinstance(model["source"], dict)
        and model["source"].get("id")
    ):
        return model["source"]["id"]


def set_reference_id(description_model, reference_id):
    if description_model.get("reference") and description_model["reference"].get("id"):
        old_reference_id = description_model["reference"]["id"]
        description_model["reference"]["id"] = reference_id
        add_msg(
            description_model["reference"],
            "info",
            {
                "code": "IUPDATEDREFERENCEID",
                "details": "Reference {} was retrieved instead of {}.".format(
                    reference_id, old_reference_id
                ),
            },
        )
    elif description_model.get("source") and description_model["source"].get("id"):
        old_reference_id = description_model["source"]["id"]
        description_model["source"]["id"] = reference_id
        add_msg(
            description_model["source"],
            "info",
            {
                "code": "IUPDATEDREFERENCEID",
                "details": "Reference {} was retrieved instead of {}.".format(
                    reference_id, old_reference_id
                ),
            },
        )
    elif description_model.get("reference"):
        description_model["reference"]["id"] = reference_id
    else:
        description_model["reference"] = {"id": reference_id}


def get_selector_id(description_model):
    """
    Get the selector ID from the description model. At the moment, no nesting
    is supported.
    :param description_model: Provided by the HGVS description parser.
    :return: The ID of the selector, if provided, otherwise None.
    """
    if (
        description_model.get("reference")
        and description_model["reference"].get("selector")
        and description_model["reference"]["selector"].get("id")
    ):
        return description_model["reference"]["selector"]["id"]
    elif (
        description_model.get("source")
        and description_model["source"].get("selector")
        and description_model["source"]["selector"].get("id")
    ):
        return description_model["source"]["selector"]["id"]


def get_coordinate_system(description_model):
    if description_model.get("coordinate_system"):
        return description_model["coordinate_system"]


def model_to_string(model):
    """
    Convert the variant description model to string.
    :param model: Dictionary holding the variant description model.
    :return: Equivalent reference string representation.
    """

    if model.get("reference"):
        reference_id = model["reference"]["id"]
    elif model.get("source"):
        reference_id = model["source"]["id"]
    selector_id = get_selector_id(model)
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


def specific_locus_to_description(specific_locus):
    """
    Convert the specific locus dictionary model to string.
    :param specific_locus: Dictionary holding the specific locus model.
    :return: Equivalent specific locus string representation.
    """
    if isinstance(specific_locus, dict):
        if specific_locus.get("id"):
            return "({})".format(specific_locus.get("id"))
    return ""


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
    deleted = inserted = ""
    if variant.get("location"):
        deleted = location_to_description(variant.get("location"))
    if variant.get("inserted"):
        inserted = inserted_to_description(variant["inserted"], sequences)
    variant_type = variant.get("type")
    if variant_type == "substitution":
        if variant.get("deleted"):
            if isinstance(variant["deleted"], dict):
                deleted += variant["deleted"]["sequence"]
            elif isinstance(variant["deleted"], list):
                deleted += variant["deleted"][0]["sequence"]
        variant_type = ">"
    elif variant_type == "deletion":
        variant_type = "del"
    elif variant_type == "deletion_insertion":
        variant_type = "delins"
    elif variant_type == "insertion":
        variant_type = "ins"
    elif variant_type == "duplication":
        variant_type = "dup"
        inserted = ""
    elif variant_type == "inversion":
        variant_type = "inv"
    elif variant_type == "equal":
        variant_type = "="
    else:
        variant_type = ""
    return "{}{}{}".format(deleted, variant_type, inserted)


def inserted_to_description(inserted, sequences):
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
        offset = "%+d" % point["offset"]["value"]
    if point.get("uncertain_offset"):
        offset = point.get("uncertain_offset")
    return "{}{}{}".format(outside_cds, position, offset)


def construct_reference(reference_id, selector_id):
    if selector_id is None:
        return {"id": reference_id}
    else:
        return {"id": reference_id, "selector": {"id": selector_id}}


def get_errors(model):
    errors = []
    if isinstance(model, list):
        for m in model:
            errors.extend(get_errors(m))
    elif isinstance(model, dict):
        if model.get("errors"):
            errors.extend(model["errors"])
        for k in model.keys():
            if k in [
                "location",
                "deleted",
                "inserted",
                "variants",
                "reference",
                "selector",
                "source",
            ]:
                errors.extend(get_errors(model[k]))
    return errors


def description_to_model(description):
    try:
        model = parse_description_to_model(description)
    except Exception as e:
        # TODO: Make it more explicit.
        model = {
            "errors": [
                {
                    "details": "Some error occured during description parsing.",
                    "raw_message": e,
                }
            ]
        }
    return model


# ---------------


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
                yield from yield_reference_ids(variant, path + [k, i])
        elif k == "inserted":
            for i, inserted in enumerate(model[k]):
                if isinstance(inserted.get("source"), dict):
                    yield inserted, tuple(path + [i])


def yield_reference_selector_ids(model, path=[]):
    for k in model.keys():
        if k in ["reference", "source"]:
            if isinstance(model[k], dict) and model[k].get("id"):
                if model[k].get("selector") and model[k]["selector"].get("id"):
                    yield model[k]["id"], model[k]["selector"]["id"], tuple(
                        path + [k, "selector", "id"]
                    )
        elif k in ["variants", "inserted"]:
            for i, sub_model in enumerate(model[k]):
                yield from yield_reference_selector_ids(sub_model, path + [k, i])


def yield_reference_selector_ids_coordinate_system(model, path=[]):
    for k in model.keys():
        if k in ["reference", "source"]:
            if isinstance(model[k], dict) and model[k].get("id"):
                if model[k].get("selector") and model[k]["selector"].get("id"):
                    yield model[k]["id"], model[k]["selector"]["id"], tuple(
                        path + [k, "selector", "id"]
                    )
        elif k in ["variants", "inserted"]:
            for i, sub_model in enumerate(model[k]):
                yield from yield_reference_selector_ids(sub_model, path + [k, i])
