def to_string(references, variants, sequences):
    return '{}{}:{}.{}'.format(
        reference_to_description(references['reference']),
        specific_locus_to_description(
            references['reference'].get('specific_locus')),
        references['reference'].get('coordinate_system'),
        variants_to_description(variants, sequences)
    )


def model_to_description(model):
    """
    Convert the variant description model to string.
    :param model: Dictionary holding the variant description model.
    :return: Equivalent reference string representation.
    """
    reference = reference_to_description(model.get('reference'))
    specific_locus = specific_locus_to_description(model.get('specific_locus'))
    if model.get('coordinate_system'):
        coordinate_system = model.get('coordinate_system') + '.'
    else:
        coordinate_system = ''
    variants = variants_to_description(model.get('variants'))
    return '{}{}:{}{}'.format(reference, specific_locus,
                              coordinate_system, variants)


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
        if specific_locus.get('type') == 'accession':
            return '({}.{})'.format(specific_locus['accession'],
                                    specific_locus['version'])
        elif specific_locus.get('type') == 'gene':
            if specific_locus.get('transcript_variant'):
                selector = '_v{}'.format(
                    specific_locus.get('transcript_variant'))
            elif specific_locus.get('protein_isoform'):
                selector = '_i{}'.format(specific_locus.get('protein_isoform'))
            else:
                selector = ''
            return '({}{})'.format((specific_locus.get('id'), selector))
        elif specific_locus.get('type') == 'lrg transcript':
            return specific_locus.get('transcript_variant')
        elif specific_locus.get('type') == 'lrg protein':
            return specific_locus.get('protein_isoform')
    return ''


def variants_to_description(variants, sequences=None):
    if isinstance(variants, list):
        variants_list = []
        for variant in variants:
            if variant.get('type') != 'equal':
                variants_list.append(variant_to_description(variant,
                                                            sequences))
        if len(variants) > 1:
            return '[{}]'.format(';'.join(variants_list))
        else:
            return variants_list[0]


def variant_to_description(variant, sequences=None):
    """
    Convert the variant dictionary model to string.
    :param p: Variant dictionary.
    :return: Equivalent variant string representation.
    """
    location = insertions = ''
    if variant.get('location'):
        location = location_to_description(variant.get('location'))
    if variant.get('insertions'):
        insertions = insertions_to_description(variant['insertions'])
    variant_type = variant.get('type')
    if variant_type == 'substitution':
        variant_type = '{}>'.format(
            sequences[variant['inserted'][0]['source']][variant['inserted'][0]['location']['position']])
    elif variant_type == 'deletion':
        variant_type = 'del'
    elif variant_type == 'equal':
        variant_type = '='
    return '{}{}{}'.format(location, variant_type, insertions)


def insertions_to_description(insertions):
    """
    Convert the insertions dictionary model to string.
    :param insertions: Insertions dictionary.
    :return: Equivalent insertions string representation.
    """
    insertions_list = []
    for i in insertions:
        if i.get('sequence'):
            insertions_list.append(i.get('sequence'))
        elif i.get('location'):
            insertions_list.append(location_to_description(i['location']))
            if i.get('inverted'):
                insertions_list[-1] += 'inv'
        elif i.get('reference_location'):
            insertions_list.append(model_to_description(i))
    if len(insertions) > 1:
        return '[{}]'.format(';'.join(insertions_list))
    else:
        return insertions_list[0]


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
        offset = '%+d' % point.get('offset')
    if point.get('uncertain_offset'):
        offset = point.get('uncertain_offset')
    return '{}{}{}'.format(outside_cds, position, offset)
