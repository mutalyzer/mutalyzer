from .util import get_start, get_end, update_position, sort_variants


def is_sorted(variants):
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


def check_start_end(variant):
    if get_start(variant['location']) <= get_end(variant['location']):
        return True
    else:
        return False


def check_substitution(variant):
    pass


def check_deletion(variant):
    pass


def check_semantics(variants):
    start_end = []
    for variant in variants:
        start_end.append(check_start_end(variant))
    output = {'sorted': is_sorted(variants),
              'overlap': is_overlap(variants),
              'start_end': start_end}

    return output


def check_semantics_reference(model, reference):
    pass
