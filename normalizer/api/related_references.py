from flask_restx import Namespace, Resource
from mutalyzer_retriever.related import get_related

from .common import errors

ns = Namespace("/")


@ns.route("/related_references/<string:reference_id>")
class RelatedReferences(Resource):
    @errors
    def get(self, reference_id):
        """Retrieve related reference ids."""
        return get_related(reference_id, timeout=10)
