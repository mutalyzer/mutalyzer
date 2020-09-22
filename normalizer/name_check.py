import copy

from extractor import describe_dna
from mutalyzer_mutator import mutate
from mutalyzer_retriever.retriever import NoReferenceError, NoReferenceRetrieved

from .checker import is_overlap, sort_variants
from .converter.to_delins import to_delins
from .converter.to_hgvs import to_hgvs_locations
from .converter.to_internal_coordinates import to_hgvs, to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .converter.variants_de_to_hgvs import de_to_hgvs
from .description import (
    description_to_model,
    get_errors,
    get_reference_id,
    get_references_from_description_model,
    model_to_string,
    variants_to_description,
    yield_reference_ids,
    yield_reference_selector_ids,
)
from .position_check import (
    check_locations,
    contains_uncertain_locations,
    identify_unsorted_locations,
)
from .protein import get_protein_description, get_protein_descriptions
from .reference import (
    get_gene_selectors,
    get_gene_selectors_hgnc,
    get_reference_id_from_model,
    get_reference_model,
    is_selector_in_reference,
)


def e_reference_not_retrieved(reference_id):
    return {
        "code": "ERETR",
        "details": "Reference {} could not be retrieved.".format(reference_id),
    }


def e_no_selector_found(reference_id, selector_id):
    return {
        "code": "ENOSELECTORFOUND",
        "details": "No {} selector found in reference {}.".format(
            selector_id, reference_id
        ),
    }


def e_selector_options(selector_id, selector_type, options):
    return {
        "code": "ESELECTOROPTIONS",
        "details": "{} selector identified as {}.".format(selector_id, selector_type),
        "options": options,
    }


def i_corrected_reference_id(original_id, corrected_id):
    return {
        "code": "ICORRECTEDREFERENCEID",
        "details": "Reference {} was retrieved instead of {}.".format(
            original_id, corrected_id
        ),
    }


def i_corrected_selector_id(original_id, corrected_id, correction_source):
    return {
        "code": "ICORRECTEDSELECTORID",
        "details": "Selector {} was corrected to {} from {}.".format(
            original_id, corrected_id, correction_source
        ),
    }


def set_by_path(dictionary, path, value):
    nested_dictionary = dictionary
    for k in path[:-1]:
        nested_dictionary = nested_dictionary[k]
    nested_dictionary[path[-1]] = value


class Description(object):
    def __init__(self, description):
        self.input_description = description
        self.input_model = description_to_model(description)
        self.corrected_model = copy.deepcopy(self.input_model)
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

        self.errors = {}
        self.infos = {}

    def _add_error(self, path, error):
        path = path
        if path in self.errors:
            self.errors[path].append(error)
        else:
            self.errors[path] = [error]

    def _add_info(self, path, info):
        path = path
        if path in self.infos:
            self.infos[path].append(info)
        else:
            self.infos[path] = [info]

    def _set_main_reference(self):
        reference_id = get_reference_id(self.corrected_model)
        if reference_id in self.references:
            self.references["reference"] = self.references[reference_id]

    def _correct_reference_id(self, path, original_id, corrected_id):
        set_by_path(self.corrected_model, path, corrected_id)
        self._add_info(path, i_corrected_reference_id(original_id, corrected_id))

    def get_references(self):
        if not self.corrected_model:
            return
        for reference_id, path in yield_reference_ids(self.input_model):
            try:
                reference_model = get_reference_model(reference_id)
            except NoReferenceError:
                self._add_error(path, e_reference_not_retrieved(reference_id))
            except NoReferenceRetrieved:
                self._add_error(path, e_reference_not_retrieved(reference_id))
            else:
                reference_id_from_model = get_reference_id_from_model(reference_model)
                if reference_id_from_model != reference_id:
                    self._correct_reference_id(
                        path, reference_id, reference_id_from_model
                    )
                    self.references[reference_id_from_model] = reference_model
                else:
                    self.references[reference_id] = reference_model
                self._set_main_reference()

    def check_selectors_in_references(self):
        if self.errors:
            return
        for reference_id, selector_id, path in yield_reference_selector_ids(
            self.corrected_model
        ):
            if not is_selector_in_reference(selector_id, self.references[reference_id]):
                self.handle_selector_not_found(reference_id, selector_id, path)

    def _correct_selector_id(self, path, original_id, corrected_id, correction_source):
        set_by_path(self.corrected_model, path, corrected_id)
        self._add_info(
            path, i_corrected_selector_id(original_id, corrected_id, correction_source)
        )

    def handle_selector_not_found(self, reference_id, selector_id, path):
        gene_selectors = get_gene_selectors(selector_id, self.references[reference_id])
        if len(gene_selectors) == 1:
            self._correct_selector_id(path, selector_id, gene_selectors[0], "gene name")
            return
        elif len(gene_selectors) > 1:
            self._add_error(
                path, e_selector_options(selector_id, "gene", gene_selectors)
            )
            return
        gene_selectors = get_gene_selectors_hgnc(
            selector_id, self.references[reference_id]
        )
        if len(gene_selectors) == 1:
            self._correct_selector_id(path, selector_id, gene_selectors[0], "gene HGNC")
            return
        elif len(gene_selectors) > 1:
            self._add_error(
                path, e_selector_options(selector_id, "gene HGNC", gene_selectors)
            )
            return
        self._add_error(path, e_no_selector_found(reference_id, selector_id))

    def check_coordinate_systems(self):
        if not self.corrected_model.get("coordinate_system"):
            self.handle_no_coordinate_system()
        if self.corrected_model.get("coordinate_system"):
            self.check_coordinate_system_consistency()

    def check_coordinate_system_consistency(self):
        pass

    def handle_no_coordinate_system(self):
        pass

    def reference_id(self):
        if self.augmented_model:
            return self.augmented_model["reference"]["id"]
        elif (
            self.input_model
            and self.input_model.get("reference")
            and self.input_model["reference"].get("id")
        ):
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
        if self.augmented_model.get("errors"):
            if (
                self.augmented_model["errors"][0].get("details")
                == "Some error occured during description parsing."
            ):
                self.augmented_model["errors"][0] = {
                    "code": "ESYNTAX",
                    "details": "A syntax error occurred.",
                }

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
            # sorted_delins_variants = sort_variants(self.delins_model["variants"])
            # print("delins variants", variants_to_description(self.delins_model["variants"]))
            # print("sorted variants:", variants_to_description(sorted_delins_variants))
            # print("are variants sorted:", sorted_delins_variants == self.delins_model["variants"])
            # print("is overlap:", is_overlap(self.delins_model["variants"]))

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
                to_selector_id=transcript_id,
            )

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
        if (
            not self.augmented_model
            or not self.internal_coordinates_model
            or get_errors(self.augmented_model)
            or get_errors(self.internal_coordinates_model)
        ):
            return
        check_locations(
            self.augmented_model,
            self.internal_coordinates_model,
            self.references[self.reference_id()],
        )

    def normalize(self):
        self.get_references()
        self.check_selectors_in_references()

        print(model_to_string(self.input_model))
        print(model_to_string(self.corrected_model))
        print(self.infos)

        # self.augment_input_model()
        # self.get_internal_coordinate_model()
        # self.check_locations()
        # self.get_internal_indexing_model()
        # identify_unsorted_locations(self.augmented_model)
        # self.get_delins_model()
        # print(identify_unsorted_locations(self.delins_model))
        # if contains_uncertain_locations(self.delins_model):
        #     return
        # self.mutate()
        # self.extract()
        # if self.de_model:
        #     self.get_de_hgvs_internal_indexing_model()
        #     self.get_de_hgvs_coordinates_model()
        #     self.get_normalized_description()
        #     self.get_equivalent_descriptions()
        #     self.get_protein_descriptions()

    def output(self):
        output = {
            "input_model": self.input_model,
        }
        if self.augmented_model:
            output["augmented_model"] = self.augmented_model
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
