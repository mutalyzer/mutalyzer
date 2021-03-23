from flask_restx import Namespace, Resource

from normalizer.viewer import view_variants

ns = Namespace("/")


@ns.route("/view_variants/<string:description>")
class ViewVariants(Resource):
    def get(self, description):
        """Visualize a variant description."""
        return view_variants(description)
