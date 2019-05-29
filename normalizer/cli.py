from .checker import is_overlap, are_sorted, check_semantics
from .validation import variants as schema_variants
from .normalizer import sanitize
from .converter import to_delins, to_hgvs, convert_indexing,\
    variants_locations_to_internal
from .to_description import to_string
import json
from hgvsparser.hgvs_parser import HgvsParser
from hgvsparser import to_model
from retriever import retriever
from crossmapper import crossmapper
from mutator.mutator import mutate
import extractor
import json
import argparse


def locations(start, end=None):
    if end:
        return {'type': 'range',
                'start': locations(start),
                'end': locations(end)}
    else:
        return {'type': 'point',
                'position': start}


def get_reference_id(description_model):
    return description_model['reference']['id']


def get_reference_model(reference_id):
    return {reference_id: retriever.retrieve(reference_id, parse=True)}


def get_reference_models(description_model):
    references_models = get_reference_model(
        description_model['reference']['id'])

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

    # reference_id = get_reference_id(description_model)
    # reference_model = retriever.retrieve(reference_id, parse=True)

    # sequences = {'reference': reference_model['sequence']}
    references = get_reference_models(description_model)

    # We need to use the crossmapper to convert the variant locations
    # to our internal indexing, using the coordinate system provided.
    # For the moment we use our trivial implementation, considering
    # only genomic coordinates as input.

    variants_internal = variants_locations_to_internal(
        variants, references, description_model['coordinate_system'])

    # variants_internal = convert_indexing(variants)

    # We need to convert the variants to delins.
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

    de_variants_sanitized = sanitize(de_variants)

    de_variants_hgvs = to_hgvs(de_variants_sanitized, sequences)

    de_variants_hgvs_indexing = convert_indexing(de_variants_hgvs, 'hgvs')

    print(to_string(description_model, de_variants_hgvs_indexing, sequences))


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('description',
                        help="HGVS variant description to be parsed")

    args = parser.parse_args()

    mutalyzer3(args.description)


if __name__ == '__main__':
    main()