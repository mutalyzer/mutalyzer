import copy

from extractor import describe_dna
from mutalyzer_mutator import mutate
from .position_check import check_points

from .converter.to_delins import to_delins
from .converter.to_hgvs import to_hgvs_locations
from .converter.to_internal_coordinates import to_internal_coordinates, to_hgvs
from .converter.to_internal_indexing import to_internal_indexing
from .converter.variants_de_to_hgvs import de_to_hgvs
from .description import (
    description_to_model,
    get_errors,
    get_references_from_description_model,
    model_to_string,
)
from .protein import get_protein_description, get_protein_descriptions


class Description(object):
    def __init__(self, description):
        self.input_description = description
        self.input_model = description_to_model(description)
        self.augmented_description = None
        self.normalized_description = None
        self.augmented_model = {}
        self.internal_coordinates_model = {}
        self.internal_indexing_model = {}
        self.delins_model = {}
        self.de_model = {}
        self.de_hgvs_internal_indexing_model = {}
        self.de_hgvs_coordinate_model = {}
        self.de_hgvs_model = {}
        self.references = {}
        self.observed_sequence = None
        self.equivalent_descriptions = None
        self.protein_descriptions = None

    def reference_id(self):
        if self.augmented_model:
            return self.augmented_model['reference']['id']
        elif (self.input_model
              and self.input_model.get("reference")
              and self.input_model["reference"].get("id")):
            return self.input_model["reference"]["id"]
        else:
            return None

    def augment_input_model(self):
        if get_errors(self.input_model):
            return
        self.augmented_model = copy.deepcopy(self.input_model)
        if get_errors(self.augmented_model):
            return
        get_references_from_description_model(self.augmented_model, self.references)
        if get_errors(self.augmented_model):
            return
        else:
            self.augmented_description = model_to_string(self.augmented_model)

    def model_parser_errors(self):
        if self.augmented_model.get('errors'):
            if (self.augmented_model["errors"][0].get("details") ==
                    "Some error occured during description parsing."):
                self.augmented_model["errors"][0] = {
                    "code": "ESYNTAX",
                    "details": "A syntax error occurred."}

    def get_internal_coordinate_model(self):
        if self.augmented_model and not get_errors(self.augmented_model):
            self.internal_coordinates_model = to_internal_coordinates(
                self.augmented_model
            )

    def get_internal_indexing_model(self):
        if self.internal_coordinates_model and not get_errors(
            self.internal_coordinates_model
        ):
            self.internal_indexing_model = to_internal_indexing(
                self.internal_coordinates_model
            )

    def get_delins_model(self):
        if self.internal_indexing_model and not get_errors(
            self.internal_indexing_model
        ):
            self.delins_model = to_delins(self.internal_indexing_model)

    def _get_sequences(self):
        """
        Retrieves a dictionary from the references with reference ids as
        keys and their corresponding sequences as values.
        """
        sequences = {k: self.references[k].sequence() for k in self.references}
        sequences["reference"] = self.references[
            self.augmented_model["reference"]["id"]
        ].sequence()
        return sequences

    def mutate(self):
        if self.delins_model and not get_errors(self.delins_model):
            self.observed_sequence = mutate(
                self._get_sequences(), self.delins_model["variants"]
            )

    def extract(self):
        if self.is_extraction_possible():
            reference_sequence = self.references[
                self.augmented_model["reference"]["id"]
            ].sequence()
            de_variants = describe_dna(reference_sequence, self.observed_sequence)
            if de_variants:
                self.de_model = {
                    "reference": copy.deepcopy(
                        self.internal_indexing_model["reference"]
                    ),
                    "coordinate_system": "i",
                    "variants": de_variants,
                }

    def is_extraction_possible(self):
        if not get_errors(self.augmented_model) and self.observed_sequence:
            return True
        return False

    def get_de_hgvs_internal_indexing_model(self):
        if self.de_model:
            self.de_hgvs_internal_indexing_model = {
                "reference": copy.deepcopy(self.internal_indexing_model["reference"]),
                "coordinate_system": "i",
                "variants": de_to_hgvs(
                    self.de_model["variants"],
                    {
                        "reference": self.references[
                            self.augmented_model["reference"]["id"]
                        ].sequence(),
                        "observed": self.observed_sequence,
                    },
                ),
            }

    def get_de_hgvs_coordinates_model(self):
        if self.augmented_model["reference"].get("selector"):
            selector_id = self.augmented_model["reference"]["selector"]["id"]
        else:
            selector_id = None
        if self.de_hgvs_internal_indexing_model:
            self.de_hgvs_model = to_hgvs_locations(
                self.de_hgvs_internal_indexing_model["variants"],
                self.references[self.augmented_model["reference"]["id"]].model,
                selector_id,
                True,
            )

    def get_normalized_description(self):
        if self.de_hgvs_model:
            self.normalized_description = model_to_string(self.de_hgvs_model)

    def get_equivalent_descriptions(self):
        if not self.de_model:
            return
        equivalent_descriptions = []

        transcript_ids = self.references[self.reference_id()].get_available_selectors()

        for transcript_id in transcript_ids[:20]:
            internal_model = to_internal_coordinates(self.de_hgvs_model)
            converted_model = to_hgvs(
                description_model=internal_model,
                to_coordinate_system=None,
                to_selector_id= transcript_id)

            equivalent_descriptions.append(model_to_string(converted_model))
        self.equivalent_descriptions = equivalent_descriptions

    def get_protein_descriptions(self):
        if self.de_model:
            references = {k: self.references[k].model for k in self.references}
            references["reference"] = self.references[self.reference_id()].model
            references["observed"] = {"sequence": {"seq": self.observed_sequence}}
            self.protein_descriptions = get_protein_descriptions(
                self.de_model["variants"], references
            )

    def check_locations(self):
        if (not self.augmented_model
            or not self.internal_coordinates_model
            or get_errors(self.augmented_model)
            or get_errors(self.internal_coordinates_model)
        ):
            return
        check_points(
            self.augmented_model,
            self.internal_coordinates_model,
            self.references[self.reference_id()])

    def normalize(self):
        self.augment_input_model()
        self.get_internal_coordinate_model()
        self.check_locations()
        self.get_internal_indexing_model()
        self.get_delins_model()
        self.mutate()
        self.extract()
        if self.de_model:
            self.get_de_hgvs_internal_indexing_model()
            self.get_de_hgvs_coordinates_model()
            self.get_normalized_description()
            self.get_equivalent_descriptions()
            self.get_protein_descriptions()

    def output(self):
        output = {
            "input_model": self.input_model,
        }
        # if self.augmented_description:
        #     output["augmented_description"] = self.augmented_description
        if self.augmented_model:
            output["augmented_model"] = self.augmented_model
        # if self.internal_coordinates_model:
        #     output["internal_coordinates_model"] = self.internal_coordinates_model
        # if self.internal_indexing_model:
        #     output["internal_indexing_model"] = self.internal_indexing_model
        # if self.normalized_description:
            output["normalized_description"] = self.normalized_description
        if self.equivalent_descriptions is not None:
            output["equivalent_descriptions"] = self.equivalent_descriptions
        if self.protein_descriptions:
            output["protein_descriptions"] = self.protein_descriptions
        return output


def normalize(description_to_normalize):
    description = Description(description_to_normalize)
    description.normalize()
    return description.output()
