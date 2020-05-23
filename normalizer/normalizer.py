from functools import lru_cache
import copy
from mutalyzer_hgvs_parser import parse_description_to_model
from retriever import retriever
from mutator.mutator import mutate
import extractor
from crossmapper import Crossmap
from .converter import de_to_hgvs, variants_locations_to_hgvs, to_delins, \
    location_to_internal, get_point_value, point_to_coordinate, coding_to_point, \
    location_to_hgvs_2, to_internal_locations
from .reference import get_selector_model, get_exon_cds_for_mrna_reference, \
    get_mol_type, get_all_selectors_exon_cds, get_transcripts_ids, \
    get_selector_model_2, get_available_selectors
from .to_description import to_string, variants_to_description
from .util import print_time_information, get_time_information
import time
import json


@lru_cache(maxsize=32)
def get_reference_model(reference_id):
    return retriever.retrieve(reference_id, parse=True)


def get_selector_id(description_model):
    """
    Get the selector ID from the description model. At the moment, no nesting
    is supported.
    :param description_model: Provided by the HGVS description parser.
    :return: The ID of the selector, if provided, otherwise None.
    """
    if description_model.get('reference') and \
            description_model['reference'].get('selector') and \
            description_model['reference']['selector'].get('id'):
        return description_model['reference']['selector']['id']


def get_coordinate_system(description_model):
    if description_model.get('coordinate_system'):
        return description_model['coordinate_system']


def string_k_v(width, key, value):
    return ' {k:<{w}} : {v}\n'.format(w=width, k=key, v=value)


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
        self._crossmap_function = None
        self._point_function = None
        self._internal_location_variants = []
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
            output += string_k_v(w, 'Reference model', 'Not retrieved')
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
            self._selector_model = get_selector_model_2(
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
            self._selector_model = get_selector_model_2(
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
                print(self._selector_id)
                if self._selector_id is None:
                    self._handle_no_description_selector()

    def _crossmapper_setup(self):
        self._mol_type = get_mol_type(
            self._reference_models[self._reference_id])

        if self._input_description_model.get('coordinate_system') is None:
            self._coordinate_system = 'g'
            self._add_warning('No coordinate system mentioned. We assumed a '
                              'genomic one.')
        else:
            self._coordinate_system = self._input_description_model.get(
                'coordinate_system')

        if self._coordinate_system == 'g':
            if self._mol_type == 'mRNA':
                self._add_error('mRNA reference for g coordinate system')
            else:
                crossmap = Crossmap()
                self._crossmap_function = crossmap.genomic_to_coordinate
                self._point_function = get_point_value
        if self._coordinate_system == 'c':
            if self._input_description_model['reference'].get('selector'):
                selector_model = get_selector_model(
                    self._reference_models[self._reference_id],
                    self._mol_type,
                    self._input_description_model['reference']['selector']['id'])
                if not selector_model:
                    selectors = get_all_selectors_exon_cds(
                        self._reference_models[self._reference_id]['model'])
                    self._add_error('Selector {} not found. '
                                    'Choose from: {}.'.format(
                        self._input_description_model['reference']['selector']['id'],
                        ', '.join(['{}'.format(i['id']) for i in selectors])))
                    return
                crossmap = Crossmap(selector_model['exons'],
                                    selector_model['cds'])
                self._crossmap_function = crossmap.coding_to_coordinate
                self._point_function = point_to_coordinate
            else:
                if 'DNA' in self._mol_type.upper():
                    selectors = get_all_selectors_exon_cds(
                        self._reference_models[self._reference_id]['model'])
                    self._add_error('No selector. Choose from: {}.'.format(
                        ', '.join(['{}'.format(i['id']) for i in selectors])))
                elif self._mol_type == 'mRNA':
                    exons, cds = get_exon_cds_for_mrna_reference(
                        self._reference_models[self._reference_id]['model'])
                    crossmap = Crossmap(exons, cds)
                    self._crossmap_function = crossmap.coding_to_coordinate
                    self._point_function = point_to_coordinate

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

    def _to_internal_locations(self):
        for variant in self._input_description_model['variants']:
            new_variant = copy.deepcopy(variant)
            new_variant['location'] = location_to_internal(
                variant['location'], variant['type'],
                self._crossmap_function, self._point_function)
            if new_variant.get('inserted'):
                for ins in new_variant['inserted']:
                    if ins.get('location'):
                        if isinstance(ins['source'], dict):
                            self._append_reference(ins['source']['id'])
                        else:
                            crossmap_function = self._crossmap_function
                            point_function = self._point_function
                        ins['location'] = location_to_internal(
                            ins['location'], None,
                            crossmap_function, point_function)

            self._internal_location_variants.append(new_variant)

    def _mutate(self):
        self._sequences = {}
        for reference_id in self._reference_models.keys():
            self._sequences[reference_id] = self._reference_models[
                reference_id]['sequence']['seq']
        self._sequences['reference'] = self._reference_models[
            self._reference_id]['sequence']['seq']
        self._sequences['observed'] = mutate(
            self._sequences, self._delins_variants)

    def fix_crossmap_coding(self, exon, coding, crossmap):
        first_exon_coding = crossmap.coordinate_to_coding(exon)
        if coding[0] == first_exon_coding[0]:
            coding = (coding[0] + coding[1], 0, coding[2])
        return coding

    def coordinate_to_coding(self, variants, selector_model):
        crossmap = Crossmap(selector_model['exons'], selector_model['cds'])
        coordinate_variants = []
        for variant in variants:
            new_variant = copy.deepcopy(variant)
            variant['location'] = location_to_hgvs_2(variant['location'],
                                                     variant['type'])
            if variant['location']['type'] == 'point':
                coding = crossmap.coordinate_to_coding(
                        variant['location']['position'])
                # compensate for the crossmaper bug
                coding = self.fix_crossmap_coding(
                    selector_model['exons'][0][0], coding, crossmap)
                coding = self.fix_crossmap_coding(
                    selector_model['exons'][-1][1], coding, crossmap)
                new_variant['location'] = coding_to_point(coding)
            if variant['location']['type'] == 'range':
                # compensate for the crossmaper bug
                first_exon_coding = crossmap.coordinate_to_coding(selector_model['exons'][0][0])
                coding_start = crossmap.coordinate_to_coding(
                        variant['location']['start']['position'])
                if coding_start[0] == first_exon_coding[0]:
                    coding_start = (coding_start[0]+coding_start[1], 0, 0)
                coding_end = crossmap.coordinate_to_coding(
                        variant['location']['end']['position'])
                if coding_end[0] == first_exon_coding[0]:
                    coding_end = (coding_end[0]+coding_end[1], 0, 0)

                new_variant['location'] = {
                    'type': 'range',
                    'start': coding_to_point(coding_start),
                    'end': coding_to_point(coding_end)
                }
            coordinate_variants.append(new_variant)
        return coordinate_variants

    def get_equivalent_descriptions(self, reference_models,
                                    coordinate_variants,
                                    overlap=None):
        equivalent_descriptions = []
        # get selectors
        transcript_ids = get_transcripts_ids(
            reference_models[self._reference_id]['model'])
        for transcript_id in transcript_ids:
            equivalent_description = self.coordinate_to_coding(
                coordinate_variants, get_selector_model(
                    reference_models[self._reference_id], self._mol_type,
                    transcript_id))
            equivalent_descriptions.append('{}({}):c.{}'.format(
                self._reference_id, transcript_id,
                variants_to_description(equivalent_description)))
        return equivalent_descriptions

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

        if self.status['errors']:
            return
        self._crossmapper_setup()
        if self.status['errors']:
            return
        # self._to_internal_locations()

        print(to_internal_locations(self._input_description_model,
                                    self._reference_models))

        if self.status['errors']:
            return
        self._delins_variants = to_delins(self._internal_location_variants)
        self._time_stamps.append(('to delins', time.perf_counter()))

        self._mutate()
        self._time_stamps.append(('mutator', time.perf_counter()))

        de_variants = extractor.describe_dna(
            self._sequences['reference'], self._sequences['observed'])

        self._time_stamps.append(('description extractor',
                                  time.perf_counter()))

        de_variants_hgvs = de_to_hgvs(de_variants, self._sequences)

        self.status['equivalent_descriptions'] = self.get_equivalent_descriptions(
            self._reference_models, de_variants_hgvs)

        de_variants_hgvs_indexing = variants_locations_to_hgvs(
            de_variants_hgvs, self._reference_models, 'g')

        if self._mol_type == 'mRNA' and self._coordinate_system == 'c':
            self._input_description_model['coordinate_system'] = 'n'
        else:
            self._input_description_model['coordinate_system'] = 'g'

        # print(variants_to_description(de_variants_hgvs))
        # print(variants_to_description(de_variants_hgvs_indexing))

        self.normalized_description = to_string(
            self._input_description_model, de_variants_hgvs_indexing,
            self._sequences)

        self.status['normalized description'] = self.normalized_description

        self._time_stamps.append(('to HGVS description', time.perf_counter()))

        self.status['time information (s)'] = get_time_information(
            self._time_stamps)


def mutalyzer3(hgvs_description):

    description = Description(hgvs_description)
    print(description)

    return {k: description.status[k] for k in description.status if description.status[k]}
