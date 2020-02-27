from normalizer.normalizer import mutalyzer3
import argparse
import json


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('description',
                        help="HGVS variant description to be parsed")

    args = parser.parse_args()

    print(json.dumps(mutalyzer3(args.description), indent=2))


if __name__ == '__main__':
    main()
