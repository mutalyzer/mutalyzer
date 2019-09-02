import copy
from .util import get_start, get_end, set_start, get_location_length, roll,\
    get_inserted_length
from crossmapper import Crossmap
import json


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


def point_to_coordinate(point):
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


def get_point_value(point):
    return point['position']


def location_to_internal(location, crossmap_function, point_function,
                         variant_type):
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
        crossmap_function=crossmap_function,
        point_function=point_function,
        variant_type=None)


def get_mol_type(reference_model):
    for part in reference_model['model']:
        if part['type'] == 'region':
            return part['qualifiers'].get('mol_type')


def get_exon_cds_for_genomic_reference(sequences, reference):
    exons = []
    cds = []
    if '_v' in reference['selector']['id']:
        gene_id = reference['selector']['id'].split('_v')[0]
        transcript_number = int(reference['selector']['id'].split('_v')[1])
        for feature in sequences[reference['id']]['model']:
            if feature['type'] == 'gene' and feature.get('sub_features') and \
                    '-' in feature['id'] and \
                    feature['id'].split('gene-')[1] == gene_id:
                rna_id = 1
                for sub_feature in feature['sub_features']:
                    if 'rna' in sub_feature['id']:
                        if rna_id == transcript_number:
                            for part in sub_feature['sub_features']:
                                if part['type'] == 'exon':
                                    exons.append(
                                        (part['start']['position'].position,
                                         part['end']['position'].position))
                                elif part['type'] == 'CDS':
                                    cds.append(
                                        part['start']['position'].position)
                                    cds.append(
                                        part['end']['position'].position)
                        rna_id += 1
    else:
        for feature in sequences[reference['id']]['model']:
            if feature['type'] == 'gene':
                if feature.get('sub_features'):
                    for sub_feature in feature['sub_features']:
                        if sub_feature['type'] == 'mRNA' and \
                                '-' in sub_feature['id'] and \
                                reference['selector']['id'] == \
                                sub_feature['id'].split('-')[1]:
                            for part in sub_feature['sub_features']:
                                if part['type'] == 'exon':
                                    exons.append(
                                        (part['start']['position'].position,
                                         part['end']['position'].position))
                                elif part['type'] == 'CDS':
                                    cds.append(
                                        part['start']['position'].position)
                                    cds.append(
                                        part['end']['position'].position)
    if len(cds) >= 2:
        cds = sorted([cds[0], cds[-1]])
    else:
        cds = []
    cds = sorted(cds)
    return sorted(exons), cds


def get_exon_cds_for_mrna_reference(sequences, reference):
    exons = []
    cds = []
    for feature in sequences[reference['id']]['model']:
        if feature['type'] == 'gene' and feature.get('sub_features'):
            for sub_feature in feature['sub_features']:
                if sub_feature['type'] == 'CDS':
                    cds.append(
                        sub_feature['start']['position'].position)
                    cds.append(
                        sub_feature['end']['position'].position)
                elif sub_feature['type'] == 'exon':
                    exons.append(
                        (sub_feature['start']['position'].position,
                         sub_feature['end']['position'].position))
    if len(cds) >= 2:
        cds = sorted([cds[0], cds[-1]])
    else:
        cds = []
    return sorted(exons), cds


def variants_locations_to_internal(variants, sequences, from_cs, reference):
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
        if not exons:
            raise Exception('No exons.')
        if not cds:
            raise Exception('No CDS.')
        crossmap = Crossmap(locations=exons, cds=cds)
        crossmap_function = crossmap.coding_to_coordinate
        point_function = point_to_coordinate
    else:
        raise(ValueError('Locations conversion from \'{}\' coordinate system '
                         'not supported.'.format(from_cs)))
    new_variants = []
    for variant in variants:
        new_variant = copy.deepcopy(variant)
        new_variant['location'] = location_to_internal(
            location=variant['location'],
            crossmap_function=crossmap_function,
            point_function=point_function,
            variant_type=variant['type'])
        if new_variant.get('inserted'):
            for insertion in new_variant['inserted']:
                if insertion.get('location'):
                    if isinstance(insertion['source'], dict):
                        if insertion.get('coordinate_system'):
                            inserted_to_internal(
                                insertion,
                                sequences,
                                insertion.get('coordinate_system'),
                                insertion.get('source'))
                    else:
                        insertion['location'] = location_to_internal(
                            location=insertion['location'],
                            crossmap_function=crossmap_function,
                            point_function=point_function,
                            variant_type=None)
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
