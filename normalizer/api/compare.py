from flask_restx import Namespace, Resource, inputs, reqparse

from normalizer.algebra import compare

ns = Namespace("/")

_args = reqparse.RequestParser()
_args.add_argument(
    "reference",
    type=str,
    help="Reference",
    default="NG_012337.3",
    required=True,
)
_args.add_argument(
    "reference_type",
    type=str,
    help="Reference",
    default="id",
    required=True,
)
_args.add_argument(
    "lhs",
    type=str,
    help="Description 1",
    default="NG_012337.3:g.100del",
    required=True,
)
_args.add_argument(
    "lhs_type",
    type=str,
    help="Description 1",
    default="hgvs",
    required=True,
)
_args.add_argument(
    "rhs",
    type=str,
    help="Description 2",
    default="NG_012337.3:g.99_101del",
    required=True,
)
_args.add_argument(
    "rhs_type",
    type=str,
    help="Description 2",
    default="hgvs",
    required=True,
)


@ns.route("/compare/")
class Compare(Resource):
    @ns.expect(_args)
    def get(self):
        """Retrieve the reference model."""
        return compare(**_args.parse_args())
