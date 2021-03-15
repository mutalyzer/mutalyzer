import logging

from flask import Blueprint, url_for
from flask_restx import Api, Resource, fields, inputs, reqparse, apidoc
from mutalyzer_hgvs_parser import parse, to_model

from normalizer.description_extractor import description_extractor
from normalizer.name_checker import name_check
from normalizer.position_converter import position_convert
from normalizer.reference import (
    get_reference_model,
    get_reference_model_segmented,
    get_selectors_ids,
)

from ..util import log_dir

logging.basicConfig(level=logging.INFO, filename=log_dir())


# Trick to make the swagger files available under "/api".
class PatchedApi(Api):
    def _register_apidoc(self, app):
        patched_api = apidoc.Apidoc(
            'restx_doc',
            'flask_restx.apidoc',
            template_folder='templates',
            static_folder='static',
            static_url_path='/api'
        )

        @patched_api.add_app_template_global
        def swagger_static(filename):
            return url_for('restx_doc.static', filename=filename)

        app.register_blueprint(patched_api)


blueprint = Blueprint("api", __name__)

api = PatchedApi(blueprint, version="1.0", title="Mutalyzer3 API")

ns = api.namespace("/")


@ns.route("/syntax_check/<string:hgvs_description>")
class SyntaxCheck(Resource):
    def get(self, hgvs_description):
        """Check the syntax correctness of a variant description."""
        try:
            parse(hgvs_description)
        except Exception:
            return "Some error occurred."
        else:
            return "Correct syntax."


@ns.route("/description_to_model/<string:hgvs_description>")
class DescriptionToModel(Resource):
    def get(self, hgvs_description):
        """Convert a variant description to its dictionary model."""
        try:
            model = to_model(hgvs_description)
        except Exception:
            model = {"errors": "Some unexpected parsing error occured."}
        return model


args_reference_model = reqparse.RequestParser()
args_reference_model.add_argument(
    "reference_id",
    type=str,
    help="Reference ID.",
    default="NG_012337.1",
    required=True,
)
args_reference_model.add_argument(
    "feature_id",
    type=str,
    help="Restrict to certain feature id.",
    default=None,
    required=False,
)
args_reference_model.add_argument(
    "siblings",
    type=inputs.boolean,
    help="Include the feature siblings.",
    default=False,
    required=False,
)
args_reference_model.add_argument(
    "ancestors",
    type=inputs.boolean,
    help="Include the feature ancestors.",
    default=True,
    required=False,
)
args_reference_model.add_argument(
    "descendants",
    type=inputs.boolean,
    help="Include the feature descendants.",
    default=True,
    required=False,
)


@ns.route("/reference_model/")
class ReferenceModel(Resource):
    @api.expect(args_reference_model)
    def get(self):
        """Retrieve the reference model."""
        args = args_reference_model.parse_args()
        return get_reference_model_segmented(**args)


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
    "position",
    type=str,
    help="Position to be converted.",
    required=False,
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
        """Generates the HGVS variant description from a reference sequence
        and an observed sequence."""
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
            selectors = get_selectors_ids(reference_model["annotations"])
            return {"reference": reference_id, "selectors": selectors}
        return {"errors": [{"code": "ERETR"}]}
