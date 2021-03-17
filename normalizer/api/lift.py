from flask_restx import Namespace, Resource, reqparse

from normalizer.lifter import lift

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


@ns.route("/lift/")
class Lift(Resource):
    @ns.expect(_args)
    def get(self):
        """Lift a description to another reference."""
        return lift(**_args.parse_args())
