import logging

from flask import Blueprint
from flask_restx import Api, Resource, fields, inputs, reqparse
from mutalyzer_hgvs_parser import parse_description, parse_description_to_model

from normalizer.description_extractor import description_extractor
from normalizer.name_checker import name_check
from normalizer.position_converter import position_convert
from normalizer.reference import get_reference_model, get_selectors_ids

from ..util import log_dir

logging.basicConfig(level=logging.INFO, filename=log_dir())

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
        except Exception:
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
        return name_check(hgvs_description)


parser = reqparse.RequestParser()
parser.add_argument(
    "description",
    type=str,
    help="Description on which the positions are considered.",
    required=False,
)
parser.add_argument(
    "reference_id",
    type=str,
    help="Reference ID on which the positions are considered.",
    required=False,
)
parser.add_argument(
    "from_selector_id",
    type=str,
    help="Selector ID from which to convert from.",
    default="",
    required=False,
)
parser.add_argument(
    "from_coordinate_system",
    type=str,
    help="Coordinate system.",
    default="",
    required=False,
)
parser.add_argument(
    "position", type=str, help="Position to be converted.", required=False,
)
parser.add_argument(
    "to_selector_id",
    type=str,
    help="Selector ID to which to convert to.",
    default="",
    required=False,
)
parser.add_argument(
    "to_coordinate_system",
    type=str,
    help="Coordinate system.",
    default="",
    required=False,
)
parser.add_argument(
    "include_overlapping",
    type=inputs.boolean,
    help="Include overlapping selectors.",
    default=False,
    required=False,
)


@ns.route("/position_convert/")
class PositionConvert(Resource):
    @api.expect(parser)
    def get(self):
        """Converts reference positions to selector orientated
        positions and vice versa."""
        args = parser.parse_args()
        return position_convert(**args)


de_parser = reqparse.RequestParser()
de_parser.add_argument(
    "reference",
    type=str,
    help="Reference sequence.",
    default="AAAATTTCCCCCGGGG",
    required=True,
)
de_parser.add_argument(
    "observed",
    type=str,
    help="Observed sequence.",
    default="AAAATTTCCCCCGGGG",
    required=True,
)


@ns.route("/description_extract/")
class DescriptionExtract(Resource):
    @api.expect(de_parser)
    def get(self):
        """Convert a position."""
        args = de_parser.parse_args()
        return description_extractor(**args)


get_selectors_model = api.model(
    "getSelectorsModel",
    {
        "reference_id": fields.String(
            description="The reference ID", required=True, example="LRG_24"
        ),
    },
)


@ns.route("/get_selectors/<string:reference_id>")
class GetSelectors(Resource):
    @api.doc(get_selectors_model)
    def get(self, reference_id):
        """Retrieve available selectors for the provided reference."""
        reference_model = get_reference_model(reference_id)
        if reference_model:
            selectors = get_selectors_ids(reference_model["model"])
            return {"reference": reference_id, "selectors": selectors}
        return {"errors": [{"code": "ERETR"}]}
