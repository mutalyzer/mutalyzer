from .util import get_start, get_end, update_position


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
