from flask_restx import Namespace, Resource, inputs, reqparse

from normalizer.mapper import map_description

from .common import errors

ns = Namespace("/")

_args = reqparse.RequestParser()

_args.add_argument(
    "description",
    type=str,
    help="Description to be mapped.",
    default="NM_003002.2:c.274G>T",
    required=True,
)

_args.add_argument(
    "reference_id",
    type=str,
    help="Reference to which the description should be mapped.",
    default="NM_003002.4",
    required=True,
)

_args.add_argument(
    "selector_id",
    type=str,
    help="Selector Id to which the description should be mapped.",
    required=False,
)

_args.add_argument(
    "slice_to",
    type=str,
    help="Slice the sequence according to the transcript or the gene locations.",
    required=False,
)

_args.add_argument(
    "filter",
    type=inputs.boolean,
    help="Filter variants that appear due to the sequences differences.",
    default=False,
    required=False,
)


@ns.route("/map/")
class Map(Resource):
    @ns.expect(_args)
    @errors
    def get(self):
        """Map a description to another reference."""
        return map_description(**_args.parse_args())
