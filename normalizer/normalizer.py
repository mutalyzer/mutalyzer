from .util import get_start, get_end, update_position, get_location_length
from hgvsparser.hgvs_parser import HgvsParser
from hgvsparser import to_model
from retriever import retriever
from mutator.mutator import mutate
import extractor
from .checker import is_overlap, are_sorted, check_semantics
from .validation import variants as schema_variants
from .converter import to_delins, de_to_hgvs,\
    variants_locations_to_internal, variants_locations_to_hgvs
from .to_description import to_string
import json
import copy
from cachetools import cached, TTLCache

# cache = TTLCache(maxsize=100, ttl=300)


def get_reference_id(description_model):
    return description_model['reference']['id']


@cached(cache={})
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
                if inserted.get('source') is not None:
                    if isinstance(inserted['source'], dict):
                        references_models.update(get_reference_model(
                            inserted['source']['id']))
    return references_models


def check_fuzzy_location():
    pass


def check_variants(variants):
    for variant in variants:
        check_fuzzy_location()


def process_reference(description_model, status):
    reference_model = retriever.retrieve(
        description_model['reference']['id'],
        parse=True)
    if reference_model is None:
        status['reference'] = 'No reference model (reference not retrieved).'


def mutalyzer3(description):
    parser = HgvsParser()
    try:
        parse_tree = parser.parse(description)
    except:
        print('fsdsdf')
        return
    description_model = to_model.convert(parse_tree)
    variants = description_model['variants']

    # print(json.dumps(description_model, indent=2))

    references = get_reference_models(description_model)

    import json
    print('{}\nvariants\n{}'.format('-' * 40, '-' * 40))
    print(json.dumps(variants, indent=2))

    variants_internal = variants_locations_to_internal(
        variants, references, description_model['coordinate_system'])

    print('{}\nvariants internal\n{}'.format('-' * 40, '-' * 40))
    print(json.dumps(variants_internal, indent=2))

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

    print('{}\nde variants\n{}'.format('-' * 40, '-' * 40))
    print(json.dumps(de_variants, indent=2))

    # print('\nde_variants:\n {}'.format(
    #     to_string(description_model, de_variants, sequences)))

    de_variants_hgvs = de_to_hgvs(de_variants, sequences)

    print('{}\nde variants hgvs\n{}'.format('-' * 40, '-' * 40))
    print(json.dumps(de_variants_hgvs, indent=2))

    # print(json.dumps(de_variants_hgvs, indent=2))

    print('\nde_variants_hgvs:\n {}'.format(
        to_string(description_model, de_variants_hgvs, sequences)))

    # de_variants_hgvs_indexing = convert_indexing(de_variants_hgvs, 'hgvs')
    de_variants_hgvs_indexing = variants_locations_to_hgvs(
        de_variants_hgvs, references, description_model['coordinate_system'])

    print('{}\nde variants hgvs indexing\n{}'.format('-' * 40, '-' * 40))
    print(json.dumps(de_variants_hgvs_indexing, indent=2))


    # print('\nde_variants_hgvs_indexing:\n {}'.format(
    #     to_string(description_model, de_variants_hgvs_indexing, sequences)))

    return to_string(description_model, de_variants_hgvs_indexing, sequences)
