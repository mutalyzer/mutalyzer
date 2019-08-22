from flask import Flask, jsonify
from normalizer.normalizer import mutalyzer3
from hgvsparser.hgvs_parser import HgvsParser
from hgvsparser import to_model

app = Flask(__name__)


@app.route("/description_to_model/<string:hgvs_description>")
def description_to_model(hgvs_description):
    """Converts a variant description to its data model."""
    parser = HgvsParser()
    parse_tree = parser.parse(hgvs_description)
    return jsonify(to_model.convert(parse_tree))


@app.route("/name_check/<string:hgvs_description>")
def name_check(hgvs_description):
    """Takes a variant description as input and checks
    whether it is correct."""
    return jsonify(mutalyzer3(hgvs_description))


app.run(port=5000, debug=True)
