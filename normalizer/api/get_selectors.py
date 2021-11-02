from flask_restx import Namespace, Resource
from mutalyzer_retriever.retriever import NoReferenceRetrieved

from normalizer.reference import retrieve_reference, get_selectors_ids

from .common import errors

ns = Namespace("/")


@ns.route("/get_selectors/<string:reference_id>")
class GetSelectors(Resource):
    @errors
    def get(self, reference_id):
        """Retrieve available selectors for the provided reference."""
        reference_model = retrieve_reference(reference_id)
        if reference_model:
            selectors = get_selectors_ids(reference_model["annotations"])
            return {"reference": reference_id, "selectors": selectors}
        return {"errors": [{"code": "ERETR"}]}
