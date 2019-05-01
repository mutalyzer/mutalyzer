import copy
from .util import get_start, get_end


def point_to_range_internal_indexing(location):
    return {'type': 'range',
            'start': {'type': 'point',
                      'position': location['position'] - 1},
            'end': {'type': 'point',
                    'position': location['position']}}


def location_to_internal_indexing(location, variant_type=None):
    new_location = copy.deepcopy(location)
    if variant_type == 'insertion':
        new_location['end']['position'] -= 1
    elif variant_type == 'deletion_insertion' and \
            new_location['type'] == 'range':
        new_location['end']['position'] -= 1
    elif location['type'] == 'range':
        new_location['start']['position'] -= 1
    elif location['type'] == 'point':
        new_location = point_to_range_internal_indexing(location)
    return new_location


def variant_to_internal_indexing(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['location'] = location_to_internal_indexing(
        variant['location'], variant['type'])
    if new_variant.get('inserted'):
        for insertion in new_variant['inserted']:
            if insertion.get('location'):
                insertion['location'] = location_to_internal_indexing(
                    insertion['location'])
    return new_variant


def to_internal_indexing(variants):
    """
    Convert variants locations to the Mutalyzer internal indexing scheme.
    """
    new_variants = []
    for variant in variants:
        new_variants.append(variant_to_internal_indexing(variant))
    return new_variants


def to_hgvs_indexing(variants):
    pass


def substitution_to_delin(variant):
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
