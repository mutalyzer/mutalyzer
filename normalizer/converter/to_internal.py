import warnings
import copy
from crossmapper import Crossmap

from normalizer.reference import get_available_selectors, get_selector_model
from normalizer.description import *
from normalizer.util import get_start, get_end, set_start, get_location_length, roll,\
    get_inserted_length, sort_location_tuples


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


def fix_selector_id(reference_models, reference_id, coordinate_system):
    available_selectors = get_available_selectors(
        reference_models[reference_id]['model'],
        coordinate_system)
    if len(available_selectors) == 0:
        raise Exception(
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
        raise Exception(
            'ENOSELECTOR: {} coordinate system used but no selector ID '
            'provided in the description. Please choose between the '
            'following selectors available in the reference: {}'.format(
                coordinate_system, available_selectors))


def get_point_value(point):
    return point['position']


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

    selector = get_selector_model(
        references[reference_id]['model'], selector_id)

    crossmap = Crossmap(
        selector['exon'], selector['cds'][0], selector['inverted'])

    return {'crossmap_function': crossmap.coding_to_coordinate,
            'point_function': point_to_x_coding}


def crossmap_noncoding_to_coordinate_setup(description, references):
    reference_id = description['reference']['id']

    selector_id = get_selector_id(description)
    if selector_id is None:
        selector_id = fix_selector_id(
            references, reference_id, 'n')

    selector = get_selector_model(
        references[reference_id]['model'], selector_id)
    selector['exon'] = sort_location_tuples(selector['exon'])
    cds = (selector['exon'][0][0], selector['exon'][-1][-1])
    crossmap = Crossmap(selector['exon'], cds, selector['inverted'])
    return {'crossmap_function': crossmap.coding_to_coordinate,
            'point_function': point_to_x_coding}



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


def location_to_internal(location, variant_type, crossmap_function,
                         point_function):
    """

    :param location:
    :param variant_type:
    :param crossmap_function:
    :param point_function:
    :return:
    """
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


def inserted_to_internal(inserted):
    crossmap = Crossmap()
    crossmap_function = crossmap.genomic_to_coordinate
    point_function = get_point_value
    return location_to_internal(
        location=inserted['location'],
        variant_type=None,
        crossmap_function=crossmap_function,
        point_function=point_function)


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
                    ins['location'] = inserted_to_internal(ins)
        new_variants.append(new_variant)

    return new_variants
