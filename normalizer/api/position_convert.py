from flask_restx import Namespace, Resource, inputs, reqparse

from normalizer.position_converter import position_convert

from .common import errors

ns = Namespace("/")

_args = reqparse.RequestParser()

_args.add_argument(
    "reference_id",
    type=str,
    help="Reference ID on which the positions are considered.",
    required=False,
)
_args.add_argument(
    "from_selector_id",
    type=str,
    help="Selector ID from which to convert from.",
    default="",
    required=False,
)
_args.add_argument(
    "from_coordinate_system",
    type=str,
    help="Coordinate system.",
    default="",
    required=False,
)
_args.add_argument(
    "position",
    type=str,
    help="Position to be converted.",
    required=False,
)
_args.add_argument(
    "to_selector_id",
    type=str,
    help="Selector ID to which to convert to.",
    default="",
    required=False,
)
_args.add_argument(
    "to_coordinate_system",
    type=str,
    help="Coordinate system.",
    default="",
    required=False,
)
_args.add_argument(
    "include_overlapping",
    type=inputs.boolean,
    help="Include overlapping selectors.",
    default=False,
    required=False,
)


@ns.route("/position_convert/")
class PositionConvert(Resource):
    @errors
    @ns.expect(_args)
    def get(self):
        """Converts reference positions to selector orientated
        positions and vice versa."""
        return position_convert(**_args.parse_args())
