from flask import Flask
from flask_restx import Resource, Api
from normalizer.normalizer import mutalyzer3, get_reference_model
from mutalyzer_hgvs_parser import parse_description, parse_description_to_model

app = Flask(__name__)
api = Api(app, version='1.0', title='Mutalyzer3 API')

ns = api.namespace('/')


@ns.route("/syntax_check/<string:hgvs_description>")
class SyntaxCheck(Resource):
    def get(self, hgvs_description):
        """Check the syntax correctness of a variant description."""
        try:
            parse_description(hgvs_description)
        except Exception:
            return 'Some error occurred.'
        else:
            return 'Correct syntax.'


@ns.route("/description_to_model/<string:hgvs_description>")
class DescriptionToModel(Resource):
    def get(self, hgvs_description):
        """Convert a variant description to its dictionary model."""
        return parse_description_to_model(hgvs_description)


@ns.route("/reference_model/<string:reference_id>")
class ReferenceModel(Resource):
    def get(self, reference_id):
        """Retrieve the reference model."""
        return get_reference_model(reference_id)


@ns.route("/name_check/<string:hgvs_description>")
class NameCheck(Resource):
    def get(self, hgvs_description):
        """Normalize a variant description."""
        return mutalyzer3(hgvs_description)
