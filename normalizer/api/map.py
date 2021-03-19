from flask_restx import Namespace, Resource, inputs, reqparse

from normalizer.mapper import map_description

ns = Namespace("/")

_args = reqparse.RequestParser()

_args.add_argument(
    "description",
    type=str,
    help="Description to be lifted.",
    default="NM_003002.2:c.274G>T",
    required=True,
)

_args.add_argument(
    "reference_id",
    type=str,
    help="Reference to which the description should be lifted.",
    default="NM_003002.4",
    required=True,
)

_args.add_argument(
    "selector_id",
    type=str,
    help="Selector Id to which the description should be lifted.",
    required=False,
)

_args.add_argument(
    "clean",
    type=inputs.boolean,
    help="Filter variants that appear due to the sequences differences.",
    default=False,
    required=False,
)


@ns.route("/map/")
class Map(Resource):
    @ns.expect(_args)
    def get(self):
        """Map a description to another reference."""
        return map_description(**_args.parse_args())
