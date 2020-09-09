import copy

from .description import (
    description_to_model,
    get_errors,
    get_references_from_description_model,
    model_to_string
)
from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .converter.to_internal import to_internal_locations
from .converter.to_delins import to_delins
from mutalyzer_mutator import mutate
from extractor import describe_dna


class Description(object):

    def __init__(self, description):
        self.description = description
        self.input_model = description_to_model(description)
        self.augmented_model = {}
        self.internal_coordinates_model = {}
        self.internal_indexing_model = {}
        self.delins_model = {}
        self.de_variants = []
        self.de_internal_indexing_model = {}
        self.de_hgvs_coordinate_model = {}
        self.references = {}
        self.observed_sequence = None

    def augment_input_model(self):
        self.augmented_model = copy.deepcopy(self.input_model)
        if get_errors(self.augmented_model):
            return
        get_references_from_description_model(self.augmented_model, self.references)
        if get_errors(self.augmented_model):
            return

    def get_internal_coordinate_model(self):
        if self.augmented_model and not get_errors(self.augmented_model):
            self.internal_coordinates_model = to_internal_coordinates(
                self.augmented_model)

    def get_internal_indexing_model(self):
        if self.internal_coordinates_model and not get_errors(self.internal_coordinates_model):
            self.internal_indexing_model = to_internal_indexing(
                self.internal_coordinates_model)

    def get_delins_model(self):
        if self.internal_indexing_model and not get_errors(self.internal_indexing_model):
            self.delins_model = to_delins(self.internal_indexing_model)

    def _get_sequences(self):
        """
        Retrieves a dictionary from the _reference_models with reference ids
        as keys and their corresponding sequences as values.
        """
        sequences = {k: self.references[k].sequence() for k in self.references}
        sequences['reference'] = self.references[self.augmented_model['reference']['id']].sequence()
        return sequences

    def _mutate(self):
        if self.delins_model and not get_errors(self.delins_model):
            self.observed_sequence = mutate(
                self._get_sequences(), self.delins_model['variants'])

    def extract(self):
        reference_sequence = self.references[self.augmented_model['reference']['id']].sequence()
        if self.observed_sequence and reference_sequence:
            self.de_variants = describe_dna(
                reference_sequence, self.observed_sequence
            )

    def get_de_internal_model(self):
        pass

    def normalize(self):
        self.augment_input_model()
        self.get_internal_coordinate_model()
        self.get_internal_indexing_model()
        self.get_delins_model()
        print(model_to_string(self.delins_model))
        self._mutate()
        self.extract()
        # print(model_to_string(self.de_variants))
        # references = {k: self.references[k].model for k in self.references.keys()}
        # old_internal = to_internal_locations(self.augmented_model, references)
        # print(model_to_string(self.internal_indexing_model))
        # print(model_to_string(old_internal))

    def output(self):
        output = {
            "input_model": self.input_model,
            "augmented_model": self.augmented_model,
            "internal_coordinates_model": self.internal_coordinates_model,
            "internal_indexing_model": self.internal_indexing_model,
            "reference_ids": list(self.references.keys())
        }
        return output


def normalize(description_to_normalize):
    description = Description(description_to_normalize)
    description.normalize()
    return description.output()
