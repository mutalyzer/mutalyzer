from flask import Flask, jsonify
from normalizer.normalizer import mutalyzer3, get_reference_model
from mutalyzer_hgvs_parser import parse_description, parse_description_to_model

app = Flask(__name__)


@app.route("/syntax_check/<string:hgvs_description>")
def syntax_check(self, hgvs_description):
    """Takes a variant description as input and checks
    whether its syntax is correct."""
    try:
        parse_description(hgvs_description)
    except Exception:
        return 'Some error occurred.'
    else:
        return 'Correct syntax.'


@app.route("/description_to_model/<string:hgvs_description>")
def description_to_model(hgvs_description):
    """Converts a variant description to its dictionary model."""
    return jsonify(parse_description_to_model(hgvs_description))


@app.route("/reference_model/<string:reference_id>")
def reference_model(self, reference_id):
    """Retrieve the reference model."""
    return get_reference_model(reference_id)


@app.route("/name_check/<string:hgvs_description>")
def name_check(hgvs_description):
    """Takes a variant description as input and checks
    whether it is correct."""
    return jsonify(mutalyzer3(hgvs_description))


app.run(port=5000, debug=True)
