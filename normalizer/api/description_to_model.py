from flask_restx import Namespace, Resource
from mutalyzer_hgvs_parser import to_model

ns = Namespace("/")


@ns.route("/description_to_model/<string:description>")
class DescriptionToModel(Resource):
    def get(self, description):
        """Convert a variant description to its dictionary model."""
        try:
            model = to_model(description)
        except Exception:
            model = {"errors": "Some unexpected parsing error occured."}
        return model
