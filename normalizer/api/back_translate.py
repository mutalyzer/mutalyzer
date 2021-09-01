from flask_restx import Namespace, Resource

from normalizer.back_translator import back_translate

ns = Namespace("/")


@ns.route("/back_translate/<string:description>")
class BackTranslate(Resource):
    def get(self, description):
        """Back translate a protein variant description."""
        return back_translate(description)
