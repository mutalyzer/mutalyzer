import argparse
import json

from normalizer.name_checker import name_check
from normalizer.normalizer import normalize


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("description", help="HGVS variant description to be parsed")

    args = parser.parse_args()

    print(normalize(args.description))

    # print(description._description_model)
    # print(description.errors)
    # print(json.dumps(mutalyzer3(args.description), indent=2))


if __name__ == "__main__":
    main()
