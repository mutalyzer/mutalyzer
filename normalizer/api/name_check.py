from flask_restx import Namespace, Resource, reqparse

from normalizer.name_checker import name_check

ns = Namespace("/")


@ns.route("/name_check/<string:description>")
class NameCheck(Resource):
    def get(self, description):
        """Normalize a variant description."""
        return name_check(description)
