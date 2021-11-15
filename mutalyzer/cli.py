import argparse
import json

from mutalyzer.name_checker import name_check

from . import usage, version


def _arg_parser():
    """
    Command line argument parsing.
    """
    parser = argparse.ArgumentParser(
        description=usage[0],
        epilog=usage[1],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("-v", action="version", version=version(parser.prog))

    parser.add_argument("description", help="HGVS variant description to be parsed")

    return parser


def main():
    parser = _arg_parser()

    try:
        args = parser.parse_args()
    except IOError as error:
        parser.error(error)

    print(json.dumps(name_check(args.description), indent=2))


if __name__ == "__main__":
    main()
