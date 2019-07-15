from .util import get_start, get_end, update_position, sort_variants


def are_sorted(variants):
    """
    Check if the provided variants list is sorted.
    """
    current_position = 0
    for variant in variants:
        if get_start(variant['location']) < current_position:
            return False
        current_position = get_start(variant['location'])
    return True


def is_overlap(variants):
    """
    Check whether there is overlap between the variants.

    TODO: Add support for fuzzy (uncertain) locations.
    """
    sorted_variants = sort_variants(variants)

    current_end_position = 0
    for variant in sorted_variants:
        if get_start(variant['location']) < current_end_position:
            return True
        else:
            current_end_position = get_end(variant['location'])
    return False


def check_location_start_end(location, indexing='internal'):
    if indexing == 'internal':
        return get_start(location) <= get_end(location)
    elif indexing == 'hgvs':
        if location['type'] == 'point':
            return True
        else:
            return get_start(location) < get_end(location)


def check_variant_start_end(variant, reference=None, indexing='internal'):
    if indexing == 'internal':
        return get_start(variant['location']) <= get_end(variant['location'])
    elif indexing == 'hgvs':
        return get_start(variant['location']) < get_end(variant['location'])


def check_substitution(variant, reference=None):
    pass


def check_deletion(variant, reference=None):
    pass


def check_semantics(variants, reference=None, indexing='internal'):
    start_end = []
    for variant in variants:
        start_end.append(check_variant_start_end(variant))
    output = {'sorted': are_sorted(variants),
              'overlap': is_overlap(variants),
              'start_end': start_end}

    return output