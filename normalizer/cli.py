from .checker import is_overlap, are_sorted, check_semantics
from .validation import variants as schema_variants
from .converter import to_delins, to_hgvs, convert_indexing
import json

from tests.commons import get_variants, SEQUENCES


def locations(start, end=None):
    if end:
        return {'type': 'range',
                'start': locations(start),
                'end': locations(end)}
    else:
        return {'type': 'point',
                'position': start}


def main():

    # print(check_semantics([variants['4_4delinsAA'], variants['10del']]))
    # print(check_semantics([variants['20_10='], variants['4_4delinsAA']]))

    variants = get_variants(['6_7delins[substitution6_7]'])

    # print(schema_variants.is_valid(variants))

    # print(check_semantics(variants, indexing='hgvs'))

    # variants_internal = convert_indexing(variants, indexing='internal')
    # print(json.dumps(variants_internal, indent=2))

    variants_delins = to_hgvs(variants, SEQUENCES)

    print(json.dumps(variants_delins, indent=2))

    # print(schema_variants.is_valid(variants_delins))
