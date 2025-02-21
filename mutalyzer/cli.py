import argparse
import json

from mutalyzer.mapper import map_description
from mutalyzer.normalizer import normalize

from . import usage, version


def _normalizer_arg_parser():
    """
    Argument parser for the normalizer.
    """
    parser = argparse.ArgumentParser(
        description=usage[0],
        epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-v", action="version", version=version(parser.prog))
    parser.add_argument("description", help="HGVS variant description to be parsed")

    return parser


def normalizer():
    """CLI interface to the normalizer."""
    parser = _normalizer_arg_parser()

    try:
        args = parser.parse_args()
    except IOError as error:
        parser.error(error)

    print(json.dumps(normalize(args.description), indent=2))


def _mapper_arg_parser():
    """
    Argument parser for the mapper.
    """
    parser = argparse.ArgumentParser(
        description=usage[0],
        epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("-v", action="version", version=version(parser.prog))
    parser.add_argument("description", help="HGVS variant description to be parsed")
    parser.add_argument(
        "--reference-id", help="reference ID onto which the description is mapped"
    )
    parser.add_argument("--selector-id", help="selector ID")
    parser.add_argument(
        "--slice-to",
        help="slice the references. Choose either 'transcript' or 'gene'.",
        choices=["transcript", "gene"],
        default=None,
    )
    parser.add_argument(
        "--filter",
        help="filter the reference sequence variant differences",
        action="store_true",
    )
    parser.add_argument(
        "--len-max", help="the maximum length of the sequences", type=int
    )
    parser.add_argument(
        "--diff-max", help="the maximum length difference of the sequences", type=int
    )

    return parser


def mapper():
    """CLI interface to the mapper."""
    parser = _mapper_arg_parser()

    try:
        args = parser.parse_args()
    except IOError as error:
        parser.error(error)

    print(json.dumps(map_description(**vars(args)), indent=2))
