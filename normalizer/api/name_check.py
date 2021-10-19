from flask_restx import Namespace, Resource, inputs, reqparse

from normalizer.name_checker import name_check

from .common import errors

ns = Namespace("/")


_args = reqparse.RequestParser()

_args.add_argument(
    "sequence",
    type=str,
    help="Reference sequence.",
    required=False,
)

_args.add_argument(
    "only_variants",
    type=inputs.boolean,
    help="The description consists only of variants.",
    default=False,
    required=False,
)


@ns.route("/name_check/<string:description>")
class NameCheck(Resource):
    @errors
    @ns.expect(_args)
    def get(self, description):
        """Normalize a variant description."""
        return name_check(description, **_args.parse_args())
