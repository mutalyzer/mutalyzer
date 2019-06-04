from .util import get_start, get_end, update_position, get_location_length
from hgvsparser.hgvs_parser import HgvsParser
from hgvsparser import to_model
from retriever import retriever
from mutator.mutator import mutate
import extractor
from .checker import is_overlap, are_sorted, check_semantics
from .validation import variants as schema_variants
from .converter import to_delins, de_to_hgvs, convert_indexing,\
    variants_locations_to_internal
from .to_description import to_string
import json
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


def get_reference_id(description_model):
    return description_model['reference']['id']


def get_reference_model(reference_id):
    return {reference_id: retriever.retrieve(reference_id, parse=True)}


def get_reference_models(description_model):
    references_models = get_reference_model(
        description_model['reference']['id'])
    if description_model['reference']['id'] is None:
        print('\n Reference not retrieved!\n')

    for variant in description_model['variants']:
        if variant.get('inserted') is not None:
            for inserted in variant.get('inserted'):
                if isinstance(inserted['source'], dict):
                    references_models.update(get_reference_model(
                        inserted['source']['id']))
    return references_models


def mutalyzer3(description):
    parser = HgvsParser()
    parse_tree = parser.parse(description)
    description_model = to_model.convert(parse_tree)
    variants = description_model['variants']

    print(json.dumps(variants, indent=2))

    references = get_reference_models(description_model)

    variants_internal = variants_locations_to_internal(
        variants, references, description_model['coordinate_system'])

    variants_delins = to_delins(variants_internal)

    sequences = {}
    for reference_id in references:
        if reference_id == description_model['reference']['id']:
            sequences['reference'] = references[reference_id]['sequence']
        else:
            sequences[reference_id] = references[reference_id]['sequence']

    observed_sequence = mutate(sequences, variants_delins)

    sequences['observed'] = observed_sequence

    de_variants = extractor.describe_dna(sequences['reference'],
                                         sequences['observed'])

    print('\nde_variants:\n {}'.format(
        to_string(description_model, de_variants, sequences)))

    de_variants_hgvs = de_to_hgvs(de_variants, sequences)

    print(json.dumps(de_variants_hgvs, indent=2))

    print('\nde_variants_hgvs:\n {}'.format(
        to_string(description_model, de_variants_hgvs, sequences)))

    de_variants_hgvs_indexing = convert_indexing(de_variants_hgvs, 'hgvs')

    print('\nde_variants_hgvs_indexing:\n {}'.format(
        to_string(description_model, de_variants_hgvs_indexing, sequences)))

    return to_string(description_model, de_variants_hgvs_indexing, sequences)
