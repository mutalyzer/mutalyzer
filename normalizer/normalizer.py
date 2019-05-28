from .util import get_start, get_end, update_position, get_location_length
import copy


def fix_start_end(variant):
    start = get_start(variant['location'])
    end = get_end(variant['location'])
    if start > end:
        update_position(variant['location'], 'start', end)
        update_position(variant['location'], 'end', start)


def fix_substitution(variant):
    pass


def normalize(variants, reference=None, observed=None):
    for variant in variants:
        fix_start_end(variant)


def sanitize_variant(variant):
    sanitized = copy.deepcopy(variant)
    if variant.get('inserted'):
        new_inserted = []
        for inserted in variant['inserted']:
            if get_location_length(inserted['location']) > 0:
                new_inserted.append(inserted)
        sanitized['inserted'] = new_inserted
    return sanitized


def sanitize(variants):
    sanitized = []
    for variant in variants:
        if variant.get('type') != 'equal':
            sanitized.append(sanitize_variant(variant))
    return sanitized