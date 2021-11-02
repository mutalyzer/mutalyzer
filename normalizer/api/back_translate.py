from flask_restx import Namespace, Resource

from normalizer.back_translator import back_translate

from .common import errors

ns = Namespace("/")


@ns.route("/back_translate/<string:description>")
class BackTranslate(Resource):
    @errors
    def get(self, description):
        """Back translate a protein variant description."""
        return back_translate(description)
