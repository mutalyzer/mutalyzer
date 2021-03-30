from .description_model import yield_sub_model
from .util import get_end, get_start, sort_variants


def are_sorted(variants):
    """
    Check if the provided variants list is sorted.
    """
    current_position = 0
    for variant in variants:
        if get_start(variant["location"]) < current_position:
            return False
        current_position = get_start(variant["location"])
    return True


def is_overlap(variants):
    """
    Check whether there is overlap between the variants.

    TODO: Add support for fuzzy (uncertain) locations.
    """
    sorted_variants = sort_variants(variants)
    positions = []
    for point, path in yield_sub_model(
        variants, ["location", "start", "end"], ["point"]
    ):
        if point.get("position") is not None:
            positions.append(point["position"])
    if positions:
        min_position = min(positions) - 1
        for variant in sorted_variants:
            if get_start(variant["location"]) <= min_position - 1:
                return True
            else:
                min_position = get_end(variant["location"])
    return False


def contains_uncertain_locations(model):
    """
    Goes through model locations to see if any is uncertain.

    :param model: description model
    :return: True when the first uncertain location if encountered
             and False if none is encountered.
    """
    for location, path in yield_sub_model(
        model, ["location", "start", "end"], ["point", "range"]
    ):
        if location.get("uncertain") or (
            location.get("offset") and location["offset"].get("uncertain")
        ):
            return True
    return False


def contains_insert_length(model):
    if model.get("variants"):
        for v in model["variants"]:
            if v.get("inserted"):
                if v["type"] != "duplication" or (
                    v["type"] == "duplication" and len(v["inserted"]) > 1
                ):
                    for i in v["inserted"]:
                        if i.get("length"):
                            return True
    return False
