from .checker import is_overlap, are_sorted, check_semantics
from .validation import variants as schema_variants
from .normalizer import sanitize
from .converter import to_delins, to_hgvs, convert_indexing
from .to_description import to_string
import json
from hgvsparser.hgvs_parser import HgvsParser
from hgvsparser.to_model import parse_tree_to_model
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
    reference = description_model['references']['reference']
    if reference.get('type') == 'genbank':
        if reference.get('version') is not None:
            return '{}.{}'.format(reference['accession'], reference['version'])
        else:
            return reference['accession']
    elif reference.get('type') == 'lrg':
        return reference['id']


def mutalyzer3(description):
    # description = 'ENSG00000157764:g.100del'

    parser = HgvsParser()
    parse_tree = parser.parse(description)
    description_model = parse_tree_to_model(parse_tree).get('model')
    variants = description_model['variants']

    reference_id = get_reference_id(description_model)
    reference_model = retriever.retrieve(reference_id, parse=True)

    sequences = {'reference': reference_model['sequence']}

    # We need to use the crossmapper to convert the variant locations
    # to our internal indexing, using the coordinate system provided.
    # For the moment we use our trivial implementation, considering
    # only genomic coordinates as input.

    variants_internal = convert_indexing(variants)

    # We need to convert the variants to delins.
    variants_delins = to_delins(variants_internal)

    sequences['observed'] = mutate(sequences, variants_delins)

    de_variants = extractor.describe_dna(sequences['reference'],
                                         sequences['observed'])

    de_variants_sanitized = sanitize(de_variants)

    de_variants_hgvs = to_hgvs(de_variants_sanitized, sequences)

    de_variants_hgvs_indexing = convert_indexing(de_variants_hgvs, 'hgvs')

    print(json.dumps(description_model['references']))
    print(to_string(description_model['references'],
                    de_variants_hgvs_indexing, sequences))


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('description',
                        help="HGVS variant description to be parsed")

    args = parser.parse_args()

    mutalyzer3(args.description)


if __name__ == '__main__':
    main()
