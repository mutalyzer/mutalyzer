import copy
from .util import get_start, get_end, set_start, get_location_length


def get_indexing_shift(indexing='internal'):
    """
    Derive the correct shift to be used when applying the provided indexing.
    """
    if indexing == 'internal':
        shift = -1
    elif indexing == 'hgvs':
        shift = +1
    return shift


def point_to_range(point_location):
    """
    Convert a point location to a range location.

    :param point_location: A point location model complying object.
    :return: A range location model complying object.
    """
    return {'type': 'range',
            'start': {'type': 'point',
                      'position': point_location['position'] - 1},
            'end': {'type': 'point',
                    'position': point_location['position']}}


def range_to_point(range_location):
    """
    Convert a point location to a range location.

    :param range_location: A point location model complying object.
    :return: A point location model complying object.
    """
    if get_start(range_location) == get_end(range_location):
        return {'type': 'point',
                'position': get_start(range_location)}
    else:
        return range_location


def convert_location_index(location, variant_type=None, indexing='internal'):
    """
    Converts a location positions according to the provided indexing.

    :param location: A location model complying object.
    :param variant_type: Required in order to be able to deal with special
    cases as, e.g., insertions.
    :param indexing: To what indexing to apply the conversion applied.
    :return: A converted location instance.
    """
    shift = get_indexing_shift(indexing)
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


def substitution_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    return new_variant


def deletion_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    return new_variant


def duplication_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    new_variant['inserted'] = [{'source': 'reference',
                               'location': copy.deepcopy(
                                   new_variant['location'])}]
    set_start(new_variant['location'], get_end(new_variant['location']))
    return new_variant


def insertion_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    return new_variant


def inversion_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    new_variant['inserted'] = [{'source': 'reference',
                                'location': copy.deepcopy(
                                    new_variant['location']),
                                'inverted': True}]
    return new_variant


def conversion_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    return new_variant


def deletion_insertion_to_delins(variant):
    return copy.deepcopy(variant)


def equal_to_delins(variant):
    """
    Only works for variants using internal indexing
    """
    new_variant = copy.deepcopy(variant)
    new_variant['inserted'] = [{'source': 'reference',
                               'location': copy.deepcopy(
                                   new_variant['location'])}]
    new_variant['type'] = 'deletion_insertion'
    return new_variant


def to_delins(variants):
    """
    Convert the variants list to its deletion insertion only equivalent. It
    considers that internal indexing is employed.
    """
    new_variants = []
    for variant in variants:
        new_variants.append(globals()[variant['type']+'_to_delins'](variant))
    return new_variants


def delins_to_substitution(variant, sequences):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'substitution'
    return new_variant


def delins_to_deletion(variant):
    return {'type': 'deletion',
            'source': 'reference',
            'location': copy.deepcopy(variant['location'])}


def delins_to_duplication(variant, sequences):
    new_variant = copy.deepcopy(variant)
    return new_variant


def delins_to_insertion(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'insertion'
    return new_variant


def delins_to_inversion(variant):
    return {'type': 'inversion',
            'source': 'reference',
            'location': copy.deepcopy(variant['location'])}


def delins_to_conversion(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    return new_variant


def delins_to_deletion_insertion(variant):
    return copy.deepcopy(variant)


def delins_to_equal(variant):
    return {'type': 'equal',
            'source': 'reference',
            'location': copy.deepcopy(variant['location'])}


def is_duplication(delins_variant, sequences):
    if get_location_length(delins_variant['location']) != 0:
        return False


def variant_to_hgvs(variant, sequences):
    if variant.get('inserted') is None:
        return delins_to_deletion(variant)

    elif len(variant.get('inserted')) == 0:
        return delins_to_deletion(variant)

    elif get_start(variant['location']) == get_end(variant['location']):
        return delins_to_insertion(variant)

    elif len(variant['inserted']) == 1:
        if variant['inserted'][0]['location'] == variant['location'] and \
                variant['inserted'][0]['source'] == 'reference':
            # Todo: what happens if the source is different than the reference
            #       but the actual sequence is the same? Should we check that?
            if variant['inserted'][0].get('inverted') is True:
                return delins_to_inversion(variant)
            return delins_to_equal(variant)
        if get_location_length(variant['location']) == \
                get_location_length(variant['inserted'][0]['location']) == 1:
            return delins_to_substitution(variant, sequences)

        if is_duplication(variant, sequences):
            return delins_to_duplication(variant)

    return delins_to_deletion_insertion(variant)


def to_hgvs(variants, sequences=None):
    """
    Convert the variants to an HGVS format (e.g., a deletion insertion of one
    nucleotide is converted to a substitution).
    """
    new_variants = []
    for variant in variants:
        new_variants.append(variant_to_hgvs(variant, sequences))
    return new_variants
