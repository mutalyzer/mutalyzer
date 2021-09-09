from flask_restx import Namespace, Resource

from mutalyzer_retriever.related import get_related

ns = Namespace("/")


@ns.route("/related_references/<string:reference_id>")
class RelatedReferences(Resource):
    def get(self, reference_id):
        """Retrieve related reference ids."""
        return get_related(reference_id, timeout=5)
