from functools import lru_cache
from mutalyzer_hgvs_parser import parse_description_to_model
from retriever import retriever
from mutator.mutator import mutate
import extractor

from .converter import to_delins
from .converter import to_internal_locations
from .converter import de_to_hgvs
from .converter import to_hgvs_locations
from .reference import get_mol_type
from .reference import get_selectors_ids
from .reference import get_selector_model
from .reference import get_available_selectors
from .description import model_to_string
from .description import get_selector_id
from .description import get_coordinate_system
from .util import get_time_information
from .util import string_k_v
import time
import json


@lru_cache(maxsize=32)
def get_reference_model(reference_id):
    return retriever.retrieve(reference_id, parse=True)


class Description(object):

    def __init__(self, description):
        self.input_description = description
        self._input_description_model = None
        self._normalized_description = None
        self._equivalent_descriptions = []
        self._reference_id = None
        self._reference_models = {}
        self._selector_id = None
        self._selector_model = None
        self._coordinate_system = None
        self._mol_type = None
        self._de_hgvs_variants = None
        self.status = {'errors': [], 'warnings': []}
        self._time_stamps = []

        self._normalize()

    def __repr__(self):
        output = '{}\n'.format(self.input_description)
        w = 20
        if self._input_description_model:
            output += string_k_v(w, 'Syntax check', 'Pass')
        else:
            output += string_k_v(w, 'Syntax check', 'Failed')
        if self._reference_id:
            output += string_k_v(w, 'Reference ID', self._reference_id)
        if self._coordinate_system:
            output += string_k_v(w, 'Coordinate system',
                                 self._coordinate_system)
        if self._reference_models:
            output += string_k_v(w, 'Reference model', 'Retrieved')
        else:
            output += string_k_v(w, 'Reference model', 'Not retrieved')
        if self._mol_type:
            output += string_k_v(w, 'Reference mol type', self._mol_type)
        else:
            output += string_k_v(w, 'Reference mol type', 'Not retrieved')
        if self._selector_id:
            output += string_k_v(w, 'Selector ID', self._selector_id)
            if self._selector_model:
                output += string_k_v(w, 'Selector model', 'Retrieved')
            else:
                output += string_k_v(w, 'Selector model', 'Not retrieved')
        else:
            output += string_k_v(w, 'Selector ID', '-')
        if self.status['errors']:
            output += ' Errors:\n'
            for error in self.status['errors']:
                output += '  - {}\n'.format(error)
        if self.status['warnings']:
            output += ' Warnings:\n'
            for warning in self.status['warnings']:
                output += '  - {}\n'.format(warning)
        return output

    def _add_error(self, error):
        self.status['errors'].append(error)

    def _add_warning(self, warning):
        self.status['warnings'].append(warning)

    def syntax_check(self):
        """
        Calls the HGVS parser and retrieves the description model.
        If successful, the self._input_description_model is populated,
        otherwise the parsing error is added to the errors list.
        """
        model = parse_description_to_model(self.input_description)
        if model.get('errors'):
            self._add_error(model['errors'])
        else:
            self._input_description_model = model

    def construct_reference(self):
        """
        Populates the instance reference attributes.
        The following steps are performed:
        - gets the reference ID from the description model.
        - retrieves the reference models.
        - retrieves the selector ID and its feature mode.
        - identifies the reference molecule type.
        """
        self._reference_id = self._input_description_model['reference']['id']
        self._append_reference(self._reference_id)
        self._selector_id = get_selector_id(self._input_description_model)
        if self._selector_id:
            self._selector_model = get_selector_model(
                self._reference_models[self._reference_id]['model'],
                self._selector_id)
        self._coordinate_system = get_coordinate_system(
            self._input_description_model)
        self._mol_type = get_mol_type(
            self._reference_models[self._reference_id])

    def _handle_no_description_selector(self):
        available_selectors = get_available_selectors(
            self._reference_models[self._reference_id]['model'],
            self._coordinate_system)
        if len(available_selectors) == 0:
            self._add_error(
                'ENOSELECTOR: {} coordinate system used but no selector ID '
                'provided in the description. In addition, there is no '
                'selector available in the reference model.'.format(
                    self._coordinate_system))
        elif len(available_selectors) == 1:
            self._add_warning(
                'WNOSELECTOR: {} coordinate system used but no selector ID '
                'provided in the description. Only {} present in the reference,'
                ' which is chosen as default.'.format(
                    self._coordinate_system, available_selectors[0]))
            self._selector_id = available_selectors[0]
            self._selector_model = get_selector_model(
                self._reference_models[self._reference_id]['model'],
                self._selector_id)
        elif len(available_selectors) > 1:
            self._add_error(
                'ENOSELECTOR: {} coordinate system used but no selector ID '
                'provided in the description. Please choose between the '
                'following selectors available in the reference: {}'.format(
                    self._coordinate_system, available_selectors))

    def check_description_reference_consistency(self):
        if self._mol_type in ['dna', 'genomic DNA']:
            if self._coordinate_system in ['c', 'n']:
                if self._selector_id is None:
                    self._handle_no_description_selector()

    def _append_reference(self, reference_id):
        """
        Appends the corresponding reference ID model to the _reference_models.
        """
        if reference_id not in self._reference_models.keys():
            reference = get_reference_model(reference_id)
            if reference is None:
                self._add_error('No reference was retrieved for {}.'.format(
                    reference_id))
            else:
                self._reference_models[reference_id] = reference

    def _mutate(self):
        self._sequences = {}
        for reference_id in self._reference_models.keys():
            self._sequences[reference_id] = self._reference_models[
                reference_id]['sequence']['seq']
        self._sequences['reference'] = self._reference_models[
            self._reference_id]['sequence']['seq']
        self._sequences['observed'] = mutate(
            self._sequences, self._delins_variants)

    def get_equivalent_descriptions(self):

        equivalent_descriptions = []

        transcript_ids = get_selectors_ids(
            self._reference_models[self._reference_id]['model'])

        for transcript_id in transcript_ids:
            equivalent_variant_model = to_hgvs_locations(
                self._de_hgvs_variants,
                self._reference_models[self._reference_id],
                transcript_id)

            equivalent_descriptions.append(model_to_string(
                equivalent_variant_model))
        return equivalent_descriptions

    def _get_normalized_description(self):
        normalized_description_model = to_hgvs_locations(
                self._de_hgvs_variants,
                self._reference_models[self._reference_id],
                self._selector_id)

        self.normalized_description = model_to_string(
            normalized_description_model)

    def _normalize(self):
        self._time_stamps.append(('initial', time.perf_counter()))

        self.syntax_check()

        if self.status['errors']:
            return

        self._time_stamps.append(('syntax parser', time.perf_counter()))

        self.construct_reference()

        self.check_description_reference_consistency()
        print(self)

        self._time_stamps.append(('retriever', time.perf_counter()))

        internal_locations_variants = to_internal_locations(
            self._input_description_model, self._reference_models)

        if self.status['errors']:
            return
        self._delins_variants = to_delins(internal_locations_variants)
        self._time_stamps.append(('to delins', time.perf_counter()))

        self._mutate()
        self._time_stamps.append(('mutator', time.perf_counter()))

        de_variants = extractor.describe_dna(
            self._sequences['reference'], self._sequences['observed'])

        self._time_stamps.append(('description extractor',
                                  time.perf_counter()))

        self._de_hgvs_variants = de_to_hgvs(de_variants, self._sequences)

        self._get_normalized_description()

        self.status['equivalent_descriptions'] = self.get_equivalent_descriptions()

        self.status['normalized description'] = self.normalized_description

        self._time_stamps.append(('to HGVS description', time.perf_counter()))

        self.status['time information (s)'] = get_time_information(
            self._time_stamps)


def mutalyzer3(hgvs_description):

    description = Description(hgvs_description)

    return {k: description.status[k] for k in description.status if description.status[k]}
