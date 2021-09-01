from flask_restx import Namespace, Resource
from mutalyzer_hgvs_parser import to_model
from mutalyzer_hgvs_parser.exceptions import UnexpectedCharacter, UnexpectedEnd
import normalizer.errors as errors

ns = Namespace("/")


@ns.route("/description_to_model/<string:description>")
class DescriptionToModel(Resource):
    def get(self, description):
        """Convert a variant description to its dictionary model."""
        try:
            model = to_model(description)
        except UnexpectedCharacter as e:
            return {"errors": [errors.syntax_uc(e)]}
        except UnexpectedEnd as e:
            return {"errors": [errors.syntax_ueof(e)]}
        return model
