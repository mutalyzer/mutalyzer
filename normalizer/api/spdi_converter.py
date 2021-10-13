from flask_restx import Namespace, Resource

from normalizer.spdi_converter import spdi_converter

ns = Namespace("/")


@ns.route("/spdi_converter/<string:description>")
class SpdiConverter(Resource):
    def get(self, description):
        """Normalize a variant description."""
        return spdi_converter(description)
