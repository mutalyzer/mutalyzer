import json
import time
from functools import lru_cache

import extractor
from mutalyzer_hgvs_parser import parse_description_to_model
from mutalyzer_mutator.mutator import mutate
from mutalyzer_retriever import retriever

from .converter import (
    de_to_hgvs,
    to_cds_coordinate,
    to_delins,
    to_hgvs_locations,
    to_internal_locations,
)
from .converter.to_internal_coordinates import location_to_internal_coordinate, to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .description import (
    get_coordinate_system,
    get_selector_id,
    model_to_string,
    variant_to_description,
    get_errors
)
from .protein import get_protein_description, get_protein_descriptions
from .reference import (
    extract_sequences,
    get_available_selectors,
    get_mol_type,
    get_selector_model,
    get_selectors_ids,
    get_reference_model
)
from .util import get_time_information, string_k_v
from .visualization import to_be_visualized
from .checker import run_checks


class NameCheck(object):
    def __init__(self, description):
        self.description = description
        self._description_model = self.to_description_model()
        self.to_internal_coordinates()
        self.errors = get_errors(self._description_model)

    def to_description_model(self):
        try:
            model = parse_description_to_model(self.description)
        except Exception as e:
            model = {"errors": [{
                "details": "Some error occured during description parsing.",
                "raw_message": e
            }]}
        return model

    def to_internal_coordinates(self):
        to_internal_coordinates(self._description_model)
