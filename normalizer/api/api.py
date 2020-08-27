from flask import Blueprint, Flask, render_template, request
from flask_restx import Api, Resource, fields, marshal, reqparse, inputs
from mutalyzer_hgvs_parser import parse_description, parse_description_to_model
from normalizer.normalizer import get_reference_model, mutalyzer3
from normalizer.position_convert import position_convert

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
        try:
            model = parse_description_to_model(hgvs_description)
        except Exception as e:
            model = {"errors": "Some unexpected parsing error occured."}
        return model


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


parser = reqparse.RequestParser()
parser.add_argument('reference_id',
                    type=str,
                    help="Reference ID.",
                    default="NG_012337.1",
                    required=True)
parser.add_argument('selector_id',
                    type=str,
                    help="Selector ID.",
                    default="NM_003002.2",
                    required=True)
parser.add_argument('position',
                    type=str,
                    help="Position to be converted.",
                    default="300",
                    required=True)
parser.add_argument('relative_to',
                    type=str,
                    help="Position relative to the reference or the selector.",
                    default="Reference",
                    required=True)
parser.add_argument('include_overlapping',
                    type=bool,
                    help="Include overlapping selectors.",
                    default=False,
                    required=False)


@ns.route("/position_convert/")
class PositionConvert(Resource):
    @api.expect(parser)
    def get(self):
        """Convert a position."""
        args = parser.parse_args()
        print(args)
        return position_convert(**args)
