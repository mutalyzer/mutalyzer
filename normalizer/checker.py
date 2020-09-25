from .description_model import yield_reference_selector_ids
from .reference import (
    get_reference_mol_type,
    get_sequence_length,
    is_selector_in_reference,
    yield_selector_ids,
)
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

    current_end_position = 0
    for variant in sorted_variants:
        if get_start(variant["location"]) <= current_end_position - 1:
            return True
        else:
            current_end_position = get_end(variant["location"])
    return False
