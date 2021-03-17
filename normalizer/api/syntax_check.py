from flask_restx import Namespace, Resource
from mutalyzer_hgvs_parser import parse

ns = Namespace("/")


@ns.route("/syntax_check/<string:description>")
class SyntaxCheck(Resource):
    def get(self, description):
        """Check the syntax correctness of a variant description."""
        try:
            parse(description)
        except Exception:
            return "Some error occurred."
        else:
            return "Correct syntax."
