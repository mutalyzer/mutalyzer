import copy
from .util import get_start, get_end, set_start, get_location_length, roll,\
    get_inserted_length
from crossmapper import Crossmap


def get_indexing_shift(indexing='internal'):
    """
    Derive the correct shift to be used when applying the provided
    indexing.
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


def location_to_internal(location, to_function, variant_type):
    new_location = copy.deepcopy(location)
    if location['type'] == 'range':
        new_location['start']['position'] = to_function(
            location['start']['position'])
        new_location['end']['position'] = to_function(
            location['end']['position'])
    elif location['type'] == 'point':
        new_location['position'] = to_function(
            location['position'])

    if variant_type == 'insertion':
        new_location['start']['position'] += 1
    elif location['type'] == 'range':
        new_location['end']['position'] += 1
    elif location['type'] == 'point':
        new_location = point_to_range(location)

    return new_location


def location_to_hgvs(location, to_function, variant_type):
    # import json
    # print('{}\nlocation\n{}'.format('-' * 40, '-' * 40))
    # print(json.dumps(location, indent=2))

    new_location = copy.deepcopy(location)
    if location['type'] == 'range':
        new_location['start']['position'] = to_function(
            location['start']['position'])
        new_location['end']['position'] = to_function(
            location['end']['position'])
    elif location['type'] == 'point':
        new_location['position'] = to_function(location['position'])

    if variant_type == 'insertion':
        new_location['start']['position'] -= 1
    elif location['type'] == 'range':
        new_location['end']['position'] -= 1
        if get_start(new_location) == get_end(new_location):
            new_location = range_to_point(new_location)

    return new_location


def insertion_to_internal():
    pass


def variants_locations_to_internal(variants, references, from_cs):
    if from_cs == 'g':
        crossmap = Crossmap()
        to_function = crossmap.genomic_to_coordinate
    else:
        raise(ValueError('Locations conversion from \'{}\' coordinate system '
                         'not supported.'.format(from_cs)))
    new_variants = []
    for variant in variants:
        new_variant = copy.deepcopy(variant)
        new_variant['location'] = location_to_internal(
            variant['location'],
            to_function,
            variant['type'])
        if new_variant.get('inserted'):
            for insertion in new_variant['inserted']:
                if insertion.get('location'):
                    if isinstance(insertion['source'], dict):
                        insertion_to_internal(insertion, references, from_cs)
                    else:
                        insertion['location'] = location_to_internal(
                            insertion['location'],
                            to_function,
                            None)
        new_variants.append(new_variant)
    return new_variants


def variants_locations_to_hgvs(variants, references, to_cs):
    if to_cs == 'g':
        crossmap = Crossmap()
        to_function = crossmap.coordinate_to_genomic
    else:
        raise ValueError('Locations conversion from \'{}\' coordinate system '
                         'not supported.'.format(to_cs))
    new_variants = []
    for variant in variants:
        new_variant = copy.deepcopy(variant)
        new_variant['location'] = location_to_hgvs(
            variant['location'],
            to_function,
            variant['type'])
        if new_variant.get('inserted'):
            for insertion in new_variant['inserted']:
                if insertion.get('location'):
                    if isinstance(insertion['source'], dict):
                        insertion_to_internal(insertion, references, to_cs)
                    else:
                        insertion['location'] = location_to_hgvs(
                            insertion['location'],
                            to_function,
                            None)
        new_variants.append(new_variant)
    return new_variants


def substitution_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    return new_variant


def deletion_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant['type'] = 'deletion_insertion'
    new_variant['inserted'] = []
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
    Convert the variants list to its deletion insertion only
    equivalent. It considers that internal indexing is employed.
    """
    new_variants = []
    for variant in variants:
        new_variants.append(globals()[variant['type']+'_to_delins'](variant))
    return new_variants


def delins_to_substitution(variant, sequences, o_index):
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


def delins_to_deletion_insertion(variant):
    return copy.deepcopy(variant)


def delins_to_equal(variant):
    return {'type': 'equal',
            'source': 'reference',
            'location': copy.deepcopy(variant['location'])}


def is_duplication(delins_variant, sequences):
    print('is duplication')


def is_deletion(delins_variant):
    if (delins_variant.get('inserted') is None) or\
            (len(delins_variant.get('inserted')) == 0):
        return True
    for inserted in delins_variant['inserted']:
        if get_location_length(inserted['location']) > 0:
            return False
    return True


def variant_to_hgvs(variant, sequences, o_index):
    if variant.get('type') != 'deletion_insertion':
        raise ValueError('The variant type is not \'deletion_insertion\', '
                         'but \'{}\'.'.format(variant.get('type')))

    if get_start(variant['location']) == get_end(variant['location']):
        return delins_to_insertion(variant)

    if (variant.get('inserted') is None) or \
            (len(variant.get('inserted')) == 0):
        return delins_to_deletion(variant)

    if len(variant['inserted']) == 1:
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
    Convert the variants to an HGVS format (e.g., a deletion insertion
    of one nucleotide is converted to a substitution).
    """
    new_variants = []
    for variant in variants:
        new_variants.append(variant_to_hgvs(variant, sequences, 0))

    return new_variants


def update_inserted_with_sequences(inserted, sequences):
    for insert in inserted:
        if insert['source'] == 'observed':
            insert['sequence'] = sequences['observed'][
                get_start(insert['location']):get_end(insert['location'])]


def de_to_hgvs(variants, sequences=None):
    """
    Convert the variants to an HGVS format (e.g., a deletion insertion
    of one nucleotide is converted to a substitution).
    """
    new_variants = []
    o_index = 0
    import json
    for variant in variants:
        if variant.get('type') == 'equal':
            o_index += get_location_length(variant['location'])
        elif variant.get('type') == 'inversion':
            o_index += get_location_length(variant['location'])
            new_variants.append(copy.deepcopy(variant))
        elif variant.get('type') == 'deletion_insertion':
            o_index += get_inserted_length(variant['inserted'])
            if get_start(variant['location']) == get_end(variant['location']):
                # delins_to_insertion
                # shift5, shift3 = roll(sequences['observed'],
                #                       o_index,
                #                       o_index + get_inserted_length(
                #                           variant['inserted']))
                # o_index += shift3
                new_variant = copy.deepcopy(variant)
                ins_length = get_location_length(
                    variant['inserted'][0]['location'])
                ins_seq = sequences['observed'][
                    get_start(variant['inserted'][0]['location']):
                    get_end(variant['inserted'][0]['location'])]
                # new_variant['location']['start']['position'] += shift3
                # new_variant['location']['end']['position'] += shift3

                if sequences['observed'][
                   get_start(variant['location']) - ins_length:
                   get_end(variant['location'])] == ins_seq:
                    new_variant['type'] = 'duplication'
                    new_variant['location']['start']['position'] = \
                        get_start(new_variant['location']) - ins_length
                else:
                    new_variant['type'] = 'insertion'

                update_inserted_with_sequences(new_variant['inserted'],
                                               sequences)
                new_variants.append(new_variant)


            elif is_deletion(variant):
                # delins_to_deletion
                new_variant = {'type': 'deletion',
                               'source': 'reference',
                               'location': copy.deepcopy(variant['location'])}
                new_variants.append(new_variant)

            elif len(variant['inserted']) == 1:
                if variant['inserted'][0]['location'] == variant['location'] \
                        and variant['inserted'][0]['source'] == 'reference':
                    if variant['inserted'][0].get('inverted') is True:
                        # delins_to_inversion
                        new_variant = {'type': 'inversion',
                                       'source': 'reference',
                                       'location': copy.deepcopy(
                                           variant['location'])}
                        new_variants.append(new_variant)
                elif get_location_length(variant['location']) == \
                        get_location_length(
                            variant['inserted'][0]['location']) == 1:
                    # delins_to_substitution
                    new_variant = copy.deepcopy(variant)
                    new_variant['type'] = 'substitution'
                    update_inserted_with_sequences(new_variant['inserted'],
                                                   sequences)
                    new_variant['deleted'] = {
                        'sequence': sequences['reference'][
                                    get_start(new_variant['location']):
                                    get_end(new_variant['location'])],
                        'source': 'reference_location'}
                    new_variants.append(new_variant)

                # elif is_duplication(variant, sequences):
                #     # delins_to_duplication
                #     new_variant = copy.deepcopy(variant)
                #     new_variants.append(new_variant)

                else:
                    # delins_to_deletion_insertion
                    new_variant = copy.deepcopy(variant)
                    update_inserted_with_sequences(new_variant['inserted'],
                                                   sequences)
                    new_variants.append(new_variant)
        else:
            raise ValueError('Unexpected variant type: \'{}\'.'.format(
                variant.get('type')))

    return new_variants
