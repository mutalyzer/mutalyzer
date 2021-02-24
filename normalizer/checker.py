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
    for point, path in yield_sub_model(variants, ["location", "start", "end"], "point"):
        if point.get("position") is not None:
            positions.append(point["position"])
    min_position = min(positions) - 1
    for variant in sorted_variants:
        if get_start(variant["location"]) <= min_position - 1:
            return True
        else:
            min_position = get_end(variant["location"])
    return False
