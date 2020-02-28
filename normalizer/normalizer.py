from functools import lru_cache
import copy
from mutalyzer_hgvs_parser import parse_description_to_model
from retriever import retriever
from mutator.mutator import mutate
import extractor
from crossmapper import Crossmap
from .converter import de_to_hgvs, variants_locations_to_hgvs, to_delins, \
    location_to_internal, get_point_value, point_to_coordinate
from .reference import get_selector_model, get_exon_cds_for_mrna_reference, \
    get_mol_type, get_all_selectors_exon_cds
from .to_description import to_string
from .util import print_time_information, get_time_information
import time


@lru_cache(maxsize=32)
def get_reference_model(reference_id):
    return retriever.retrieve(reference_id, parse=True)


class Description(object):

    def __init__(self, description):
        self.status = {'errors': [], 'warnings': []}
        self.description = description
        self._description_model = None
        self._reference_id = None
        self._reference_models = {}
        self._mol_type = None
        self._coordinate_system = None
        self._crossmap_function = None
        self._point_function = None
        self._time_stamps = []

        self._normalize()

    def _add_error(self, error):
        self.status['errors'].append(error)

    def _add_warning(self, warning):
        self.status['warnings'].append(warning)

    def _parse(self):
        model = parse_description_to_model(self.description)
        if model.get('errors'):
            self._add_error(model['errors'])
        else:
            self._description_model = model

    def _crossmapper_setup(self):

        self._mol_type = get_mol_type(
            self._reference_models[self._reference_id])

        if self._description_model.get('coordinate_system') is None:
            self._coordinate_system = 'g'
            self._add_warning('No coordinate system mentioned. We assumed a '
                              'genomic one.')
        else:
            self._coordinate_system = self._description_model.get(
                'coordinate_system')

        if self._coordinate_system == 'g':
            if self._mol_type == 'mRNA':
                self._add_error('mRNA reference for g coordinate system')
            else:
                crossmap = Crossmap()
                self._crossmap_function = crossmap.genomic_to_coordinate
                self._point_function = get_point_value
        if self._coordinate_system == 'c':
            if self._description_model['reference'].get('selector'):
                selector_model = get_selector_model(
                    self._reference_models[self._reference_id],
                    self._mol_type,
                    self._description_model['reference']['selector']['id'])
                if not selector_model:
                    selectors = get_all_selectors_exon_cds(
                        self._reference_models[self._reference_id]['model'])
                    self._add_error('Selector {} not found. '
                                    'Choose from: {}.'.format(
                        self._description_model['reference']['selector']['id'],
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
        valid = True
        if reference_id not in self._reference_models.keys():
            reference = get_reference_model(reference_id)
            if reference is None:
                self._add_error('No reference was retrieved for {}.'.format(
                    reference_id))
                return False
            if reference['model'] is None:
                self._add_error('No model for {}.'.format(reference_id))
                valid = False
            if reference['sequence'] is None:
                self._add_error('No sequence for {}.'.format(reference_id))
                valid = False
            if valid:
                self._reference_models[reference_id] = reference
        return valid

    def _to_internal_locations(self):
        self._internal_location_variants = []
        for variant in self._description_model['variants']:
            new_variant = copy.deepcopy(variant)
            new_variant['location'] = location_to_internal(
                variant['location'], self._crossmap_function,
                self._point_function, variant['type'])
            if new_variant.get('inserted'):
                for ins in new_variant['inserted']:
                    if ins.get('location'):
                        if isinstance(ins['source'], dict):
                            if not self._append_reference(ins['source']['id']):
                                break
                        else:
                            crossmap_function = self._crossmap_function
                            point_function = self._point_function
                        ins['location'] = location_to_internal(
                            ins['location'], crossmap_function,
                            point_function, None)

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

    def _normalize(self):
        self._time_stamps.append(('initial', time.perf_counter()))

        self._parse()
        if self.status['errors']:
            return

        self._time_stamps.append(('syntax parser', time.perf_counter()))

        self._reference_id = self._description_model['reference']['id']
        self._append_reference(self._reference_id)

        self._time_stamps.append(('retriever', time.perf_counter()))

        if self.status['errors']:
            return
        self._crossmapper_setup()
        if self.status['errors']:
            return
        self._to_internal_locations()
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

        de_variants_hgvs_indexing = variants_locations_to_hgvs(
            de_variants_hgvs, self._reference_models, 'g')

        if self._mol_type == 'mRNA' and self._coordinate_system == 'c':
            self._description_model['coordinate_system'] = 'n'
        else:
            self._description_model['coordinate_system'] = 'g'

        self.normalized_description = to_string(
            self._description_model, de_variants_hgvs_indexing,
            self._sequences)

        self.status['normalized_description'] = self.normalized_description

        self._time_stamps.append(('to HGVS description', time.perf_counter()))

        self.status['time_information (s)'] = get_time_information(
            self._time_stamps)


def mutalyzer3(hgvs_description):

    description = Description(hgvs_description)

    return {k: description.status[k] for k in description.status if description.status[k]}
