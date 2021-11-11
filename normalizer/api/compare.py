from flask_restx import Namespace, Resource, reqparse

from normalizer.algebra import compare

from .common import errors

ns = Namespace("/")

_args = reqparse.RequestParser()
_args.add_argument(
    "reference",
    type=str,
    help="Reference id or sequence",
    required=False,
)
_args.add_argument(
    "reference_type",
    type=str,
    help="Reference type",
    required=False,
)
_args.add_argument(
    "lhs",
    type=str,
    help="Left hand side operator",
    required=True,
)
_args.add_argument(
    "lhs_type",
    type=str,
    help="Left hand side operator type",
    required=True,
)
_args.add_argument(
    "rhs",
    type=str,
    help="Right hand side operator",
    required=True,
)
_args.add_argument(
    "rhs_type",
    type=str,
    help="Right hand side operator type",
    required=True,
)


@ns.route("/compare/")
class Compare(Resource):
    @ns.expect(_args)
    @errors
    def get(self):
        """Compute the relation between variants."""
        return compare(**_args.parse_args())
