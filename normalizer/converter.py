import copy
from .util import get_start, get_end


def get_index_shift(indexing='internal'):
    """
    Derive the correct shift to be used when applying the provided indexing.
    """
    if indexing == 'internal':
        shift = -1
    elif indexing == 'hgvs':
        shift = +1
    return shift


def point_to_range(point):
    """
    Convert a point location to a range location.

    :param point: A point location model complying object.
    :return: A range location model complying object.
    """
    return {'type': 'range',
            'start': {'type': 'point',
                      'position': point['position'] - 1},
            'end': {'type': 'point',
                    'position': point['position']}}


def range_to_point(range):
    """
    Convert a point location to a range location.

    :param range: A point location model complying object.
    :return: A point location model complying object.
    """
    if get_start(range) == get_end(range):
        return {'type': 'point',
                'position': get_start(range)}
    else:
        return range


def convert_location_index(location, variant_type=None, indexing='internal'):
    """
    Converts a location positions according to the provided indexing.

    :param location: A location model complying object.
    :param variant_type: Required in order to be able to deal with special
    cases as, e.g., insertions.
    :param indexing: To what indexing to apply the conversion applied.
    :return: A converted location instance.
    """
    shift = get_index_shift(indexing)
    new_location = copy.deepcopy(location)
    if variant_type == 'insertion':
        new_location['end']['position'] += shift
    elif location['type'] == 'range':
        new_location['start']['position'] += shift
        if indexing == 'hgvs':
            new_location = range_to_point(new_location)
    elif location['type'] == 'point':
        new_location = point_to_range(location)
    return new_location


def convert_variant_indexing(variant, indexing='internal'):
    """
    Convert a variant locations to the internal/HGVS indexing scheme.
    """
    new_variant = copy.deepcopy(variant)
    new_variant['location'] = convert_location_index(
        location=variant['location'],
        variant_type=variant['type'],
        indexing=indexing)
    if new_variant.get('inserted'):
        for insertion in new_variant['inserted']:
            if insertion.get('location'):
                insertion['location'] = convert_location_index(
                    location=insertion['location'],
                    indexing=indexing)
    return new_variant


def convert_indexing(variants, indexing='internal'):
    """
    Converts variants locations to the internal/HGVS indexing scheme.

    :param variants: List with variants to be converted.
    :param indexing: To what indexing to convert to.
    :return: Converted variants list instance.
    """
    new_variants = []
    for variant in variants:
        new_variants.append(convert_variant_indexing(variant, indexing))
    return new_variants


def substitution_to_delin(variant):
    """
    Only works for
    :param variant:
    :return:
    """
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    return new_variant


def deletion_to_delin(variant):
    new_variant = []
    return new_variant


def duplication_to_delin(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    return new_variant


def insertion_to_delin():
    new_variant = []
    return new_variant


def to_delins(variants):
    """
    Convert to only deletion insertion variants.
    """
    new_variants = []
    for variant in variants:
        if variant['type'] == 'deletion_insertion':
            new_variants.append(variant)
        elif variant['type'] == 'substitution':
            new_variants.append(substitution_to_delin(variant))
        elif variant['type'] == 'deletion':
            new_variants.append((deletion_to_delin(variant)))
        elif variant['type'] == 'duplication':
            new_variants.append((duplication_to_delin(variant)))
    return new_variants


def to_hgvs(variants):
    """
    Convert the variants to an HGVS format (e.g., a deletion insertion of one
    nucleotide is converted to a substitution).
    """
    pass
