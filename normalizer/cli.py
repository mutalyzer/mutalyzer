from normalizer.normalizer import mutalyzer3
import json
import argparse


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('description',
                        help="HGVS variant description to be parsed")

    args = parser.parse_args()

    mutalyzer3(args.description)


if __name__ == '__main__':
    main()
