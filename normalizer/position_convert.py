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
    to_hgvs,
    to_internal_locations,
    to_internal
)
from .description import location_to_description
from .protein import get_protein_description, get_protein_descriptions
from .reference import (
    extract_sequences,
    get_available_selectors,
    get_mol_type,
    get_selector_model,
    get_selectors_ids,
)
from .util import get_time_information, string_k_v
from .visualization import to_be_visualized
from .checker import run_checks

from mutalyzer_crossmapper import Genomic, NonCoding, Coding


def crossmap_to_x_setup(selector_model, mol_type, relative_to):
    if mol_type in ['genomic DNA', 'dna']:
        if relative_to == "Reference":
            crossmap = Genomic()
            return {
                "crossmap_function": crossmap.genomic_to_coordinate,
                "point_function": to_internal.get_point_value,
            }
        elif relative_to == "Selector" and selector_model["type"] in ["mRNA"]:
            crossmap = Coding(selector_model["exon"], selector_model["cds"][0],
                              selector_model["inverted"])
            return {
                "crossmap_function": crossmap.coding_to_coordinate,
                "point_function": to_internal.point_to_x_coding,
            }
        elif relative_to == "selector" and selector_model["type"] in ["ncRNA"]:
            crossmap = Coding(selector_model["exon"], selector_model["cds"][0],
                              selector_model["inverted"])
            return {
                "crossmap_function": crossmap.coding_to_coordinate,
                "point_function": to_internal.point_to_x_coding,
            }


def crossmap_to_hgvs_setup(selector_model, mol_type, relative_to):
    if relative_to == "Reference":
        if selector_model['type'] in ["mRNA"]:
            crossmap = Coding(selector_model["exon"], selector_model["cds"][0],
                              selector_model["inverted"])
            return {
                "crossmap_function": crossmap.coordinate_to_coding,
                "point_function": to_hgvs.coding_to_point,
                "degenerate": True
            }


def position_convert(reference_id, selector_id, position, relative_to):
    reference_model = retriever.retrieve(reference_id, parse=True)
    if reference_model:
        mol_type = get_mol_type(reference_model)
    else:
        return 'Reference not retrieved'
    location_model = parse_description_to_model(position, start_rule='location')
    if reference_model:
        selector_model = get_selector_model(reference_model["model"], selector_id)
        crossmap = crossmap_to_x_setup(selector_model, mol_type, relative_to)
        internal = to_internal.point_to_coding(location_model, **crossmap)
        crossmap = crossmap_to_hgvs_setup(selector_model, mol_type, relative_to)
        hgvs = to_hgvs.point_to_hgvs(internal, **crossmap)
        return {
            "reference": {
                "id": reference_id,
                "position": position
            },
            "selector": {
                "id": selector_id,
                "position": location_to_description(hgvs)
            }
        }



