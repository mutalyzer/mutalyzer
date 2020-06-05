

def to_string(references, variants, sequences):
    """
    Convert the variant description model to string.
    :param references:
    :param variants:
    :param sequences:
    :return: Equivalent reference string representation.
    """
    return '{}{}:{}.{}'.format(
        references['reference']['id'],
        specific_locus_to_description(
            references['reference'].get('selector')),
        references['coordinate_system'],
        variants_to_description(variants, sequences))


def reference_to_description(reference):
    """
    Convert the reference dictionary model to string.
    :param reference: Dictionary holding the reference model.
    :return: Equivalent reference string representation.
    """
    version = ''
    if isinstance(reference, dict):
        if reference.get('type') == 'genbank':
            accession = reference.get('accession')
            if reference.get('version'):
                version = '.{}'.format(reference['version'])
        elif reference.get('type') == 'lrg':
            accession = reference.get('id')
    return '{}{}'.format(accession, version)


def specific_locus_to_description(specific_locus):
    """
    Convert the specific locus dictionary model to string.
    :param specific_locus: Dictionary holding the specific locus model.
    :return: Equivalent specific locus string representation.
    """
    if isinstance(specific_locus, dict):
        if specific_locus.get('id'):
            return '({})'.format(specific_locus.get('id'))
    return ''


def variants_to_description(variants, sequences=None):
    if isinstance(variants, list):
        variants_list = []
        for variant in variants:
            variants_list.append(variant_to_description(variant, sequences))
        if len(variants_list) > 1:
            return '[{}]'.format(';'.join(variants_list))
        elif len(variants_list) == 1:
            return variants_list[0]


def variant_to_description(variant, sequences=None):
    """
    Convert the variant dictionary model to string.
    :return: Equivalent variant string representation.
    """
    deleted = inserted = ''
    if variant.get('location'):
        deleted = location_to_description(variant.get('location'))
    if variant.get('inserted'):
        inserted = inserted_to_description(variant['inserted'], sequences)
    variant_type = variant.get('type')
    if variant_type == 'substitution':
        if variant.get('deleted'):
            deleted += variant['deleted']['sequence']
        variant_type = '>'
    elif variant_type == 'deletion':
        variant_type = 'del'
    elif variant_type == 'deletion_insertion':
        variant_type = 'delins'
    elif variant_type == 'insertion':
        variant_type = 'ins'
    elif variant_type == 'duplication':
        variant_type = 'dup'
        inserted = ''
    elif variant_type == 'inversion':
        variant_type = 'inv'
    elif variant_type == 'equal':
        variant_type = '='
    return '{}{}{}'.format(deleted, variant_type, inserted)


def inserted_to_description(inserted, sequences):
    """
    Convert the insertions dictionary model to string.
    :param inserted: Insertions dictionary.
    :return: Equivalent insertions string representation.
    """
    descriptions = []
    for insert in inserted:
        if insert.get('sequence'):
            descriptions.append(insert['sequence'])
        elif insert.get('location'):
            descriptions.append(location_to_description(insert['location']))
            if insert.get('inverted'):
                descriptions[-1] += 'inv'
        elif insert.get('reference_location'):
            descriptions.append(model_to_string(insert))
    if len(inserted) > 1:
        return '[{}]'.format(';'.join(descriptions))
    else:
        return descriptions[0]


def location_to_description(location):
    """
    Convert the location dictionary model to string.
    :param location: Location dictionary.
    :return: Equivalent location string representation.
    """
    if location['type'] == 'point':
        return point_to_description(location)
    if location['type'] == 'range':
        if location.get('uncertain'):
            return '({}_{})'.format(
                point_to_description(location.get('start')),
                point_to_description(location.get('end')))
        else:
            start = location_to_description(location.get('start'))
            end = location_to_description(location.get('end'))
            return '{}_{}'.format(start, end)


def point_to_description(point):
    """
    Convert the position dictionary model to string.
    :param point: Position dictionary.
    :return: Equivalent position string representation.
    """
    outside_cds = offset = ''
    if point.get('outside_cds'):
        if point['outside_cds'] == 'downstream':
            outside_cds = '*'
        elif point['outside_cds'] == 'upstream':
            outside_cds = '-'
    if point.get('uncertain'):
        position = '?'
    else:
        position = str(point.get('position'))
    if point.get('offset'):
        offset = '%+d' % point['offset']['value']
    if point.get('uncertain_offset'):
        offset = point.get('uncertain_offset')
    return '{}{}{}'.format(outside_cds, position, offset)
