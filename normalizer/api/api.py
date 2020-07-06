from mutalyzer_crossmapper import Crossmap
from flask import Blueprint, Flask, render_template, request
from flask_restx import Api, Resource, fields, marshal
from mutalyzer_hgvs_parser import parse_description, parse_description_to_model
from normalizer.normalizer import get_reference_model, mutalyzer3

blueprint = Blueprint("api", __name__)

api = Api(blueprint, version="1.0", title="Mutalyzer3 API")

ns = api.namespace("/")


@ns.route("/syntax_check/<string:hgvs_description>")
class SyntaxCheck(Resource):
    def get(self, hgvs_description):
        """Check the syntax correctness of a variant description."""
        try:
            parse_description(hgvs_description)
        except Exception:
            return "Some error occurred."
        else:
            return "Correct syntax."


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


crossmap_input = api.model(
    "CrossmapInput",
    {
        "exons": fields.List(fields.Integer, description="Exons"),
        "cds": fields.List(fields.Integer, description="CDS"),
    },
)


@ns.route("/crossmap/")
class Crossmap(Resource):
    @api.doc(body=crossmap_input)
    def post(self):
        print(request.data)
        print(marshal(request.data, crossmap_input))
        return "OK"
