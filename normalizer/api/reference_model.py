from flask_restx import Namespace, Resource, inputs, reqparse

from normalizer.reference import get_reference_model_segmented

from .common import errors

ns = Namespace("/")

_args = reqparse.RequestParser()
_args.add_argument(
    "reference_id",
    type=str,
    help="Reference ID.",
    default="NG_012337.1",
    required=True,
)
_args.add_argument(
    "feature_id",
    type=str,
    help="Restrict to certain feature id.",
    default=None,
    required=False,
)
_args.add_argument(
    "siblings",
    type=inputs.boolean,
    help="Include the feature siblings.",
    default=False,
    required=False,
)
_args.add_argument(
    "ancestors",
    type=inputs.boolean,
    help="Include the feature ancestors.",
    default=True,
    required=False,
)
_args.add_argument(
    "descendants",
    type=inputs.boolean,
    help="Include the feature descendants.",
    default=True,
    required=False,
)


@ns.route("/reference_model/")
class ReferenceModel(Resource):
    @ns.expect(_args)
    @errors
    def get(self):
        """Retrieve the reference model."""
        return get_reference_model_segmented(**_args.parse_args())
