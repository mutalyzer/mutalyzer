import json
import copy
import warnings
from .util import get_start, get_end, set_start, get_location_length, roll,\
    get_inserted_length, sort_location_tuples
from crossmapper import Crossmap
from .reference import *
from .description import *


def point_to_range(point_location):
    """
    Convert a point location to a range location.

    :param point_location: A point location model complying object.
    :return: A range location model complying object.
    """
    return {'type': 'range',
            'start': {'type': 'point',
                      'position': point_location['position']},
            'end': {'type': 'point',
                    'position': point_location['position'] + 1}}


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


def get_noncoding_downstream_position(crossmap, selector_model):
    pass


def point_to_coordinate(point, selector_model):
    position = point['position']
    if point.get('outside_cds'):
        if point['outside_cds'] == 'upstream':
            section = 0
            position = -1 * position
        elif point['outside_cds'] == 'downstream':
            section = 2
    else:
        section = 1
    if point.get('offset'):
        offset = point['offset']['value']
    else:
        offset = 0
    return position, offset, section


def point_to_x_coding(point):
    position = point['position']
    if point.get('outside_cds'):
        if point['outside_cds'] == 'upstream':
            section = 0
            position = -1 * position
        elif point['outside_cds'] == 'downstream':
            section = 2
    else:
        section = 1
    if point.get('offset'):
        offset = point['offset']['value']
    else:
        offset = 0
    return position, offset, section


def point_to_x_noncoding(point):
    position = point['position']
    if point.get('offset'):
        offset = point['offset']['value']
    else:
        offset = 0
    if point.get('outside_cds'):
        if point['outside_cds'] == 'upstream':
            position = 1
            offset = -1 * position
        elif point['outside_cds'] == 'downstream':
            offset = position
            # TODO: Include the orientation.
    return position, offset


def x_coding_to_point(coding):
    position, offset, section = coding
    point = {'type': 'point',
             'position': position}

    if section == 0:
        point['outside_cds'] = 'upstream'
        point['position'] *= -1
    elif section == 2:
        point['outside_cds'] = 'downstream'

    if offset != 0:
        point['offset'] = {'value': offset}

    return point


def coding_to_point(coding):
    position, offset, section = coding
    point = {'type': 'point',
             'position': position}

    if section == 0:
        point['outside_cds'] = 'upstream'
        point['position'] *= -1
    elif section == 2:
        point['outside_cds'] = 'downstream'

    if offset != 0:
        point['offset'] = {'value': offset}

    return point


def x_noncoding_to_point(noncoding):
    position, offset = noncoding
    point = {'type': 'point',
             'position': position}

    if offset != 0:
        point['offset'] = {'value': offset}

    return point


def get_point_value(point):
    return point['position']


def location_to_internal(location, variant_type,
                         crossmap_function, point_function):
    new_location = copy.deepcopy(location)
    if location['type'] == 'range':
        new_location['start']['position'] = crossmap_function(
            point_function(location['start']))
        new_location['end']['position'] = crossmap_function(
            point_function(location['end']))
    elif location['type'] == 'point':
        new_location['position'] = crossmap_function(
            point_function(location))
    if variant_type == 'insertion':
        new_location['start']['position'] += 1
    elif location['type'] == 'range':
        new_location['end']['position'] += 1
    elif location['type'] == 'point':
        new_location = point_to_range(new_location)
    return new_location


def location_to_hgvs(location, to_function, variant_type):
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


def location_to_hgvs_2(location, variant_type):
    new_location = copy.deepcopy(location)

    if variant_type == 'insertion':
        new_location['start']['position'] -= 1
    elif location['type'] == 'range':
        new_location['end']['position'] -= 1
        if get_start(new_location) == get_end(new_location):
            new_location = range_to_point(new_location)

    return new_location


def inserted_to_internal(inserted, sequences, from_cs, reference):
    if from_cs == 'g':
        crossmap = Crossmap()
        crossmap_function = crossmap.genomic_to_coordinate
        point_function = get_point_value
    elif from_cs == 'c':
        mol_type = get_mol_type(sequences[reference['id']])
        if mol_type == 'genomic DNA':
            exons, cds = get_exon_cds_for_genomic_reference(sequences,
                                                            reference)
        if mol_type == 'mRNA':
            exons, cds = get_exon_cds_for_mrna_reference(sequences, reference)
        crossmap = Crossmap(locations=exons, cds=cds)
        crossmap_function = crossmap.coding_to_coordinate
        point_function = point_to_coordinate
    else:
        raise (ValueError('Locations conversion from \'{}\' coordinate system '
                          'not supported.'.format(from_cs)))
    inserted['location'] = location_to_internal(
        location=inserted['location'],
        variant_type=None,
        crossmap_function=crossmap_function,
        point_function=point_function)


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
            location=variant['location'],
            to_function=to_function,
            variant_type=variant['type'])
        if new_variant.get('inserted'):
            for insertion in new_variant['inserted']:
                if insertion.get('location'):
                    if isinstance(insertion['source'], dict):
                        inserted_to_internal(insertion, references, to_cs)
                    else:
                        insertion['location'] = location_to_hgvs(
                            location=insertion['location'],
                            to_function=to_function,
                            variant_type=None)
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
    for variant in variants:
        if variant.get('type') == 'equal':
            o_index += get_location_length(variant['location'])
        elif variant.get('type') == 'inversion':
            o_index += get_location_length(variant['location'])
            new_variants.append(copy.deepcopy(variant))
        elif variant.get('type') == 'deletion_insertion':
            if get_start(variant['location']) == get_end(variant['location']):
                # delins_to_insertion
                shift5, shift3 = roll(sequences['observed'],
                                      o_index + 1,
                                      o_index + get_inserted_length(
                                          variant['inserted']))
                ins_seq = sequences['observed'][
                    o_index:o_index + get_inserted_length(variant['inserted'])]
                o_index += shift3
                new_variant = copy.deepcopy(variant)
                ins_length = get_location_length(
                    variant['inserted'][0]['location'])
                new_variant['location']['start']['position'] += shift3
                new_variant['location']['end']['position'] += shift3

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
                shift5, shift3 = roll(
                    sequences['reference'],
                    new_variant['location']['start']['position'] + 1,
                    new_variant['location']['end']['position'])
                new_variant['location']['start']['position'] += shift3
                new_variant['location']['end']['position'] += shift3
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
                else:
                    # delins_to_deletion_insertion
                    new_variant = copy.deepcopy(variant)
                    update_inserted_with_sequences(new_variant['inserted'],
                                                   sequences)
                    new_variants.append(new_variant)
            o_index += get_inserted_length(variant['inserted'])

        else:
            raise ValueError('Unexpected variant type: \'{}\'.'.format(
                variant.get('type')))

    return new_variants


def fix_selector_id(reference_models, reference_id, coordinate_system):
    available_selectors = get_available_selectors(
        reference_models[reference_id]['model'],
        coordinate_system)
    if len(available_selectors) == 0:
        raise(
            'ENOSELECTOR: {} coordinate system used but no selector ID '
            'provided in the description. In addition, there is no '
            'selector available in the reference model.'.format(
                coordinate_system))
    elif len(available_selectors) == 1:
        warnings.warn(
            'WNOSELECTOR: {} coordinate system used but no selector ID '
            'provided in the description. Only {} present in the reference,'
            ' which is chosen as default.'.format(
                coordinate_system, available_selectors[0]))
        return available_selectors[0]
    elif len(available_selectors) > 1:
        raise(
            'ENOSELECTOR: {} coordinate system used but no selector ID '
            'provided in the description. Please choose between the '
            'following selectors available in the reference: {}'.format(
                coordinate_system, available_selectors))


def crossmap_genomic_to_coordinate_setup():
    crossmap = Crossmap()
    return {'crossmap_function': crossmap.genomic_to_coordinate,
            'point_function': get_point_value}


def crossmap_coding_to_coordinate_setup(description, references):
    reference_id = description['reference']['id']

    selector_id = get_selector_id(description)
    if selector_id is None:
        selector_id = fix_selector_id(
            references, reference_id, 'c')

    selector = get_selector_model_2(
        references[reference_id]['model'], selector_id)

    crossmap = Crossmap(selector['exon'], selector['cds'], selector['inverted'])
    return {'crossmap_function': crossmap.coding_to_coordinate,
            'point_function': point_to_x_coding}


def crossmap_noncoding_to_coordinate_setup(description, references):
    reference_id = description['reference']['id']

    selector_id = get_selector_id(description)
    if selector_id is None:
        selector_id = fix_selector_id(
            references, reference_id, 'n')

    selector = get_selector_model_2(
        references[reference_id]['model'], selector_id)
    selector['exon'] = sort_location_tuples(selector['exon'])
    crossmap = Crossmap(selector['exon'], None, selector['inverted'])
    return {'crossmap_function': crossmap.noncoding_to_coordinate,
            'point_function': point_to_x_noncoding}


def crossmap_to_x_setup(description, references):
    """
    Returns a crossmap instance able to convert from the coordinate system
    provided in the description model to the to internal system (crossmap
    coordinate).
    :param description: Description model.
    :param references: References models.
    :return: {'crossmap_function' : ..., 'point_function': ...}
    """
    coordinate_system = get_coordinate_system(description)
    if coordinate_system is None:
        warnings.warn('No coordinate system, we assume g.')
        # TODO: Improve message
        coordinate_system = 'g'
        # TODO: Update also the description model.

    if coordinate_system == 'g':
        crossmap = crossmap_genomic_to_coordinate_setup()
    elif coordinate_system == 'c':
        crossmap = crossmap_coding_to_coordinate_setup(description, references)
    elif coordinate_system == 'n':
        crossmap = crossmap_noncoding_to_coordinate_setup(description,
                                                          references)
    else:
        raise Exception('Unsupported coordinate system: {}.'.format(
            coordinate_system))

    return crossmap


def to_internal_locations(description, references):
    """
    Converts the variant locations present in the description model to the
    internal coordinate system.

    :param description: Description model, dictionary.
    :param references: References models, dictionary.
    :return: Variants with locations in the internal coordinate system.
    """

    crossmap = crossmap_to_x_setup(description, references)
    new_variants = []

    for variant in description['variants']:
        new_variant = copy.deepcopy(variant)
        new_variant['location'] = location_to_internal(
            variant['location'], variant['type'], **crossmap)
        if new_variant.get('inserted'):
            for ins in new_variant['inserted']:
                if ins.get('location'):
                    pass
        new_variants.append(new_variant)

    return new_variants