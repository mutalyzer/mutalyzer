import logging

from flask import Blueprint, url_for
from flask_restx import Api, apidoc

from ..util import log_dir
from .back_translate import ns as ns_back_translate
from .compare import ns as ns_compare
from .description_extract import ns as ns_description_extract
from .description_to_model import ns as ns_description_to_mdel
from .get_selectors import ns as ns_get_selectors
from .map import ns as ns_map
from .mutate import ns as ns_mutate
from .name_check import ns as ns_name_check
from .position_convert import ns as ns_position_convert
from .reference_model import ns as ns_reference_model
from .related_references import ns as ns_related_references
from .spdi_converter import ns as ns_spdi_converter
from .view_variants import ns as ns_view_variants

logging.basicConfig(level=logging.INFO, filename=log_dir())


# Trick to make the swagger files available under "/api".
class PatchedApi(Api):
    def _register_apidoc(self, app):
        patched_api = apidoc.Apidoc(
            "restx_doc",
            "flask_restx.apidoc",
            template_folder="templates",
            static_folder="static",
            static_url_path="/api",
        )

        @patched_api.add_app_template_global
        def swagger_static(filename):
            return url_for("restx_doc.static", filename=filename)

        app.register_blueprint(patched_api)


blueprint = Blueprint("api", __name__)

api = PatchedApi(blueprint, version="1.0", title="Mutalyzer3 API")

api.add_namespace(ns_compare)
api.add_namespace(ns_map)
api.add_namespace(ns_mutate)
api.add_namespace(ns_name_check)
api.add_namespace(ns_description_to_mdel)
api.add_namespace(ns_reference_model)
api.add_namespace(ns_related_references)
api.add_namespace(ns_position_convert)
api.add_namespace(ns_description_extract)
api.add_namespace(ns_get_selectors)
api.add_namespace(ns_view_variants)
api.add_namespace(ns_spdi_converter)
api.add_namespace(ns_back_translate)
