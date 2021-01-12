import copy

from extractor import describe_dna
from mutalyzer_hgvs_parser import parse_description_to_model
from mutalyzer_hgvs_parser.exceptions import UnexpectedCharacter, UnexpectedEnd
from mutalyzer_mutator import mutate
from mutalyzer_retriever.retriever import NoReferenceError, NoReferenceRetrieved

from .checker import is_overlap, sort_variants
from .converter.to_delins import to_delins
from .converter.to_hgvs_coordinates import to_hgvs_locations
from .converter.to_hgvs_indexing import to_hgvs_indexing
from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .converter.variants_de_to_hgvs import de_to_hgvs
from .description_model import (
    get_locations_start_end,
    get_reference_id,
    location_to_description,
    model_to_string,
    point_to_description,
    variant_to_description,
    variants_to_description,
    yield_point_locations_for_main_reference,
    yield_reference_ids,
    yield_reference_selector_ids,
    yield_reference_selector_ids_coordinate_system,
)
from .position_check import (
    check_locations,
    contains_uncertain_locations,
    identify_unsorted_locations,
)
from .protein import get_protein_description, get_protein_descriptions
from .reference import (
    get_coordinate_system_from_reference,
    get_coordinate_system_from_selector_id,
    get_gene_selectors,
    get_gene_selectors_hgnc,
    get_only_selector_id,
    get_protein_selector_model,
    get_reference_id_from_model,
    get_reference_model,
    get_selectors_ids,
    get_sequence_length,
    is_only_one_selector,
    is_selector_in_reference,
    update_start_end,
    yield_overlap_ids,
    yield_selector_ids,
)
from .util import get_end, get_start, set_by_path


def e_mismatch(input_description, model_description):
    return {
        "code": "EMISMATCH",
        "details": "Model description {} different than the input description {}.".format(
            model_description, input_description
        ),
    }


def e_reference_not_retrieved(reference_id, path):
    return {
        "code": "ERETR",
        "details": "Reference {} could not be retrieved.".format(reference_id),
        "paths": [path],
    }


def e_no_selector_found(reference_id, selector_id, path):
    return {
        "code": "ENOSELECTORFOUND",
        "details": "No {} selector found in reference {}.".format(
            selector_id, reference_id
        ),
        "paths": [path],
    }


def e_selector_options(selector_id, selector_type, options, path):
    return {
        "code": "ESELECTOROPTIONS",
        "details": "{} selector identified as {}.".format(selector_id, selector_type),
        "options": options,
        "paths": [path],
    }


def e_no_coordinate_system(path):
    return {
        "code": "ENOCOORDINATESYSTEM",
        "details": "A coordinate system is required.",
        "paths": [path],
    }


def e_coordinate_system_mismatch(
    coordinate_system, mismatch_id, mismatch_coordinate_system, path
):
    return {
        "code": "ECOORDINATESYSTEMMISMATCH",
        "details": "Coordinate system {} does not match with {} "
        "{} coordinate system. ".format(
            coordinate_system, mismatch_id, mismatch_coordinate_system
        ),
        "paths": [path],
    }


def e_out_of_boundary_lesser(position, path):
    return {
        "code": "EOUTOFBOUNDARY",
        "details": "Position {} is lesser than 1.".format(
            point_to_description(position)
        ),
        "paths": [path],
    }


def e_out_of_boundary_greater(point, sequence_length, path):
    return {
        "code": "EOUTOFBOUNDARY",
        "details": "Position {} is greater than the sequence {} length.".format(
            point_to_description(point), sequence_length
        ),
        "paths": [path],
    }


def e_insertion_range_not_consecutive(location, path):
    return {
        "code": "EINSERTIONRANGE",
        "details": "Range positions {} not consecutive in insertion location.".format(
            location_to_description(location)
        ),
        "paths": [path],
    }


def e_insertion_location_not_range(point, path):
    return {
        "code": "EINSERTIONLOCATION",
        "details": "Insertion location {} is not range.".format(
            point_to_description(point)
        ),
        "paths": [path],
    }


def e_repeat_reference_sequence_length(path):
    return {
        "code": "EREPEATREFERENCELENGTH",
        "details": "Reference sequence length not a multiple of the inserted sequence length.",
        "paths": [path],
    }


def e_repeat_sequences_mismatch(reference_sequence, repeat_sequence, path):
    return {
        "code": "EREPEATMISMATCH",
        "details": "Reference sequence {} does not contain the repeat sequence {}.".format(
            reference_sequence, repeat_sequence
        ),
        "paths": [path],
    }


def e_repeat_not_supported(variant, path):
    return {
        "code": "EREPEATUNSUPPORTED",
        "details": "Repeat variant {} not supported.".format(
            variant_to_description(variant)
        ),
        "paths": [path],
    }


def i_corrected_reference_id(original_id, corrected_id, path):
    return {
        "code": "ICORRECTEDREFERENCEID",
        "details": "Reference {} was retrieved instead of {}.".format(
            corrected_id, original_id
        ),
        "paths": [path],
    }


def i_corrected_selector_id(original_id, corrected_id, correction_source, path):
    return {
        "code": "ICORRECTEDSELECTORID",
        "details": "Selector {} was corrected to {} from {}.".format(
            original_id, corrected_id, correction_source
        ),
        "paths": [path],
    }


def i_corrected_coordinate_system(coordinate_system, correction_source, path):
    return {
        "code": "ICORRECTEDCOORDINATESYSTEM",
        "details": "Coordinate system corrected to {} from {}.".format(
            coordinate_system, correction_source
        ),
        "paths": [path],
    }


def check_errors(fn):
    def wrapper(self):
        if not self.errors:
            fn(self)
        if self.errors and self.stop_on_errors:
            raise Exception(str(self.errors))

    return wrapper


class Description(object):
    def __init__(self, description=None, description_model=None, stop_on_error=False):

        self.errors = []
        self.infos = []

        self.input_description = description if description else None
        self.input_model = description_model if description_model else {}
        self._check_input()

        self.stop_on_errors = stop_on_error

        self.corrected_model = copy.deepcopy(self.input_model)
        self.internal_coordinates_model = {}
        self.internal_indexing_model = {}
        self.delins_model = {}
        self.de_model = {}
        self.de_hgvs_internal_indexing_model = {}
        self.de_hgvs_coordinate_model = {}
        self.de_hgvs_model = {}
        self.normalized_description = None

        self.references = {}

        self.equivalent_descriptions = None
        self.protein_descriptions = None

    def _check_input(self):
        if self.input_description and not self.input_model:
            self._convert_description_to_model()
        elif self.input_description is None and self.input_model:
            # TODO: check the input_model
            self.input_description = model_to_string(self.input_model)
        elif self.input_description and self.input_model:
            # TODO: check the input_model
            model_description = model_to_string(self.input_model)
            if self.input_description != model_description:
                e_mismatch(model_description, self.input_description)

    def _add_error(self, error):
        self.errors.append(error)

    def _add_info(self, info):
        self.infos.append(info)

    @check_errors
    def _convert_description_to_model(self):
        try:
            self.input_model = parse_description_to_model(self.input_description)
        except UnexpectedCharacter as e:
            self._add_error(
                dict(
                    {"code": "ESYNTAXUC", "details": "Unexpected character."},
                    **e.serialize()
                )
            )
        except UnexpectedEnd as e:
            self._add_error(
                dict(
                    {"code": "ESYNTAXUEOF", "details": "Unexpected end of input."},
                    **e.serialize()
                )
            )

    def _set_main_reference(self):
        reference_id = get_reference_id(self.corrected_model)
        if reference_id in self.references:
            self.references["reference"] = self.references[reference_id]

    def _correct_reference_id(self, path, original_id, corrected_id):
        set_by_path(self.corrected_model, path, corrected_id)
        self._add_info(i_corrected_reference_id(original_id, corrected_id, path))

    @check_errors
    def retrieve_references(self):
        """
        Populate the references
        :return:
        """
        if not self.corrected_model:
            return
        for reference_id, path in yield_reference_ids(self.input_model):
            try:
                reference_model = get_reference_model(reference_id)
                # print(get_reference_model.cache_info())
            except NoReferenceError:
                self._add_error(e_reference_not_retrieved(reference_id, [path]))
            except NoReferenceRetrieved:
                self._add_error(e_reference_not_retrieved(reference_id, [path]))
            else:
                reference_id_in_model = get_reference_id_from_model(reference_model)
                if reference_id_in_model != reference_id:
                    self._correct_reference_id(
                        path, reference_id, reference_id_in_model
                    )
                    self.references[reference_id_in_model] = reference_model
                else:
                    self.references[reference_id] = reference_model
                self._set_main_reference()

    @check_errors
    def _check_selectors_in_references(self):
        for reference_id, selector_id, path in yield_reference_selector_ids(
            self.corrected_model
        ):
            if not is_selector_in_reference(selector_id, self.references[reference_id]):
                self.handle_selector_not_found(reference_id, selector_id, path)

    def _correct_selector_id(self, path, original_id, corrected_id, correction_source):
        set_by_path(self.corrected_model, path, corrected_id)
        self._add_info(
            i_corrected_selector_id(original_id, corrected_id, correction_source, path)
        )

    def handle_selector_not_found(self, reference_id, selector_id, path):
        gene_selectors = get_gene_selectors(selector_id, self.references[reference_id])
        if len(gene_selectors) == 1:
            self._correct_selector_id(path, selector_id, gene_selectors[0], "gene name")
            return
        elif len(gene_selectors) > 1:
            self._add_error(
                e_selector_options(selector_id, "gene", gene_selectors, path)
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
                e_selector_options(selector_id, "gene HGNC", gene_selectors, path)
            )
            return
        self._add_error(e_no_selector_found(reference_id, selector_id, path))

    @check_errors
    def _check_coordinate_systems(self):
        for (
            c_s,
            c_s_path,
            r_id,
            r_path,
            s_id,
            s_path,
        ) in yield_reference_selector_ids_coordinate_system(
            copy.deepcopy(self.corrected_model)
        ):
            if c_s is None:
                self._handle_no_coordinate_system(c_s_path, r_id, s_id)

    def _correct_coordinate_system(self, coordinate_system, path, correction_source):
        set_by_path(self.corrected_model, path, coordinate_system)
        self._add_info(
            i_corrected_coordinate_system(coordinate_system, correction_source, path)
        )

    def _handle_no_coordinate_system(self, c_s_path, r_id, s_id):
        if s_id:
            c_s = get_coordinate_system_from_selector_id(self.references[r_id], s_id)
            if c_s:
                self._correct_coordinate_system(c_s, c_s_path, s_id + " selector")
                return
        c_s = get_coordinate_system_from_reference(self.references[r_id])
        if c_s:
            self._correct_coordinate_system(c_s, c_s_path, r_id + " reference")
            return
        self._add_error(e_no_coordinate_system(c_s_path))

    def _correct_selector_id_from_coordinate_system(self, r_id_path, selector_id):
        path = tuple(list(r_id_path[:-1]) + ["selector"])
        set_by_path(self.corrected_model, path, {"id": selector_id})
        self._add_info(
            i_corrected_selector_id("", selector_id, "coordinate system", path)
        )

    @check_errors
    def _check_coordinate_system_consistency(self):
        for (
            c_s,
            c_s_path,
            r_id,
            r_path,
            s_id,
            s_path,
        ) in yield_reference_selector_ids_coordinate_system(
            copy.deepcopy(self.corrected_model)
        ):
            if s_id:
                s_c_s = get_coordinate_system_from_selector_id(
                    self.references[r_id], s_id
                )
                if s_c_s == c_s:
                    return
                else:
                    self._add_error(
                        e_coordinate_system_mismatch(c_s, s_id, s_c_s, c_s_path)
                    )
                    return
            r_c_s = get_coordinate_system_from_reference(self.references[r_id])
            if (r_c_s == c_s) and (c_s in ["g", "m"]):
                return
            else:
                if is_only_one_selector(self.references[r_id], c_s):
                    self._correct_selector_id_from_coordinate_system(
                        r_path, get_only_selector_id(self.references[r_id], c_s)
                    )
                    return
                else:
                    self._add_error(
                        e_coordinate_system_mismatch(c_s, r_id, r_c_s, c_s_path)
                    )

    @check_errors
    def _construct_internal_coordinate_model(self):
        self.internal_coordinates_model = to_internal_coordinates(
            self.corrected_model, self.references
        )

    @check_errors
    def _construct_internal_indexing_model(self):
        self.internal_indexing_model = to_internal_indexing(
            self.internal_coordinates_model
        )

    @check_errors
    def _construct_delins_model(self):
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
        sequences = {k: self.references[k]["sequence"]["seq"] for k in self.references}
        return sequences

    @check_errors
    def _mutate(self):
        if self.delins_model:
            observed_sequence = mutate(
                self._get_sequences(), self.delins_model["variants"]
            )
            self.references["observed"] = {"sequence": {"seq": observed_sequence}}

    @check_errors
    def _extract(self):
        de_variants = describe_dna(
            self.references["reference"]["sequence"]["seq"],
            self.references["observed"]["sequence"]["seq"],
        )
        if de_variants:
            self.de_model = {
                "reference": copy.deepcopy(self.internal_indexing_model["reference"]),
                "coordinate_system": "i",
                "variants": de_variants,
            }

    @check_errors
    def _construct_de_hgvs_internal_indexing_model(self):
        if self.de_model:
            self.de_hgvs_internal_indexing_model = {
                "reference": copy.deepcopy(self.internal_indexing_model["reference"]),
                "coordinate_system": "i",
                "variants": de_to_hgvs(
                    self.de_model["variants"],
                    self._get_sequences(),
                ),
            }

    @check_errors
    def _construct_de_hgvs_coordinates_model(self):
        if self.corrected_model["reference"].get("selector"):
            selector_id = self.corrected_model["reference"]["selector"]["id"]
        else:
            selector_id = None
        if self.de_hgvs_internal_indexing_model:
            hgvs_indexing = to_hgvs_indexing(self.de_hgvs_internal_indexing_model)
            self.de_hgvs_model = to_hgvs_locations(
                hgvs_indexing,
                self.references,
                self.corrected_model["coordinate_system"],
                selector_id,
                True,
            )

    def _construct_normalized_description(self):
        if self.de_hgvs_model:
            self.normalized_description = model_to_string(self.de_hgvs_model)

    def _construct_equivalent_descriptions(self):
        if not self.de_model:
            return
        equivalent_descriptions = {}
        get_locations_start_end(self.de_hgvs_internal_indexing_model)

        start_limit, end_limit = update_start_end(
            self.references["reference"],
            *get_locations_start_end(self.de_hgvs_internal_indexing_model)
        )

        internal_model = to_internal_coordinates(self.de_hgvs_model, self.references)

        if (
            get_coordinate_system_from_reference(self.references["reference"])
            == "g"
            != self.corrected_model["coordinate_system"]
        ):
            equivalent_descriptions["g"] = [
                model_to_string(
                    to_hgvs_locations(
                        internal_model=internal_model,
                        references=self.references,
                        to_coordinate_system="g",
                        to_selector_id=None,
                        degenerate=True,
                    )
                )
            ]

        for selector in yield_overlap_ids(
            self.references["reference"], start_limit, end_limit
        ):
            converted_model = to_hgvs_locations(
                internal_model=internal_model,
                references=self.references,
                to_coordinate_system=None,
                to_selector_id=selector["id"],
                degenerate=True,
            )
            c_s = converted_model["coordinate_system"]
            if not equivalent_descriptions.get(c_s):
                equivalent_descriptions[c_s] = []

            if converted_model["coordinate_system"] == "c":
                protein_selector_model = get_protein_selector_model(
                    self.references["reference"]["annotations"], selector["id"]
                )
                if protein_selector_model:
                    equivalent_descriptions[c_s].append(
                        (
                            model_to_string(converted_model),
                            get_protein_description(
                                self.de_model["variants"],
                                self.references,
                                protein_selector_model,
                            ),
                        )
                    )
                else:
                    equivalent_descriptions[c_s].append(
                        model_to_string(converted_model)
                    )
            else:
                equivalent_descriptions[c_s].append(model_to_string(converted_model))

        self.equivalent_descriptions = equivalent_descriptions

    def _construct_protein_descriptions(self):
        if self.de_model:
            self.protein_descriptions = get_protein_descriptions(
                self.de_model["variants"], self.references
            )

    def check_locations(self):
        pass

    def _check_location_boundaries(self, path):
        v = self.internal_coordinates_model["variants"][path[1]]
        v_r = self.input_model["variants"][path[1]]
        if v["location"]["type"] == "point" and not v.get("uncertain"):
            if (
                get_sequence_length(self.references, "reference")
                < v["location"]["position"]
            ):
                self._add_error(
                    e_out_of_boundary_greater(
                        v_r["location"],
                        get_sequence_length(self.references, "reference"),
                        path,
                    )
                )
            elif v["location"]["position"] < 0:
                self._add_error(e_out_of_boundary_lesser(v_r["location"], path))
        if v["location"]["type"] == "range":
            if v["location"]["start"]["type"] == "point" and not v.get("uncertain"):
                if (
                    get_sequence_length(self.references, "reference")
                    < v["location"]["start"]["position"]
                ):
                    self._add_error(
                        e_out_of_boundary_greater(
                            v_r["location"]["start"],
                            get_sequence_length(self.references, "reference"),
                            path + ["start"],
                        )
                    )
                elif v["location"]["start"]["position"] < 0:
                    self._add_error(
                        e_out_of_boundary_lesser(
                            v_r["location"]["start"], path + ["start"]
                        )
                    )
            if v["location"]["end"]["type"] == "point" and not v.get("uncertain"):
                if (
                    get_sequence_length(self.references, "reference")
                    < v["location"]["end"]["position"]
                ):
                    self._add_error(
                        e_out_of_boundary_greater(
                            v_r["location"]["end"],
                            get_sequence_length(self.references, "reference"),
                            path + ["end"],
                        )
                    )
                elif v["location"]["end"]["position"] < 0:
                    self._add_error(
                        e_out_of_boundary_lesser(v_r["location"]["end"], path + ["end"])
                    )

    def _check_location_range(self, path):
        v = self.internal_coordinates_model["variants"][path[1]]
        v_r = self.input_model["variants"][path[1]]
        if v["location"]["type"] == "point" and not v.get("uncertain"):
            if (
                get_sequence_length(self.references, "reference")
                < v["location"]["position"]
            ):
                self._add_error(
                    e_out_of_boundary_greater(
                        v_r["location"],
                        get_sequence_length(self.references, "reference"),
                        path,
                    )
                )
            elif v["location"]["position"] < 0:
                self._add_error(e_out_of_boundary_lesser(v_r["location"], path))

    def _check_insertion_location(self, path):
        v = self.internal_coordinates_model["variants"][path[1]]
        v_r = self.input_model["variants"][path[1]]
        if v["location"]["type"] == "point" and not v["location"].get("uncertain"):
            self._add_error(e_insertion_location_not_range(v_r["location"], path))

        if v["location"]["type"] == "range" and not v["location"].get("uncertain"):
            if (
                abs(
                    v["location"]["start"]["position"]
                    - v["location"]["end"]["position"]
                )
                != 1
            ):
                self._add_error(
                    e_insertion_range_not_consecutive(v_r["location"], path)
                )

    def _check_repeat(self, path):
        v = self.input_model["variants"][path[1]]
        v_i = self.internal_indexing_model["variants"][path[1]]
        if v.get("inserted") and len(v.get("inserted")) == 1:
            inserted = v["inserted"][0]
            if inserted.get("sequence") and inserted.get("source") == "description":
                repeat_seq = inserted["sequence"]
            # TODO: get the sequence from a reference slice
            else:
                self._add_error(e_repeat_not_supported(v, path))
                return

            ref_seq = self.references["reference"]["sequence"]["seq"][
                get_start(v_i) : get_end(v_i)
            ]

            if len(ref_seq) % len(repeat_seq) != 0:
                self._add_error(e_repeat_reference_sequence_length(path))
            elif (len(ref_seq) // len(repeat_seq)) * repeat_seq != ref_seq:
                self._add_error(e_repeat_sequences_mismatch(ref_seq, repeat_seq, path))
        else:
            # TODO: Convert to delins and switch to warning?
            self._add_error(e_repeat_not_supported(v, path))

    @check_errors
    def check(self):
        for i, variant in enumerate(self.internal_coordinates_model["variants"]):
            if variant.get("location"):
                path = ["variants", i, "location"]
                self._check_location_boundaries(path)

            if variant.get("type") == "insertion":
                self._check_insertion_location(["variants", i])

            if variant.get("type") == "repeat":
                self._check_repeat(["variants", i])

    def to_internal_coordinate_model(self):
        self.retrieve_references()

        self._check_selectors_in_references()
        self._check_coordinate_systems()
        self._check_coordinate_system_consistency()

        self._construct_internal_coordinate_model()

    def normalize(self):
        self.retrieve_references()

        self._check_selectors_in_references()
        self._check_coordinate_systems()
        self._check_coordinate_system_consistency()

        self._construct_internal_coordinate_model()
        self._construct_internal_indexing_model()
        self._construct_delins_model()

        if contains_uncertain_locations(self.delins_model):
            return

        self.check()

        self._mutate()
        self._extract()
        if self.de_model:
            self._construct_de_hgvs_internal_indexing_model()
            self._construct_de_hgvs_coordinates_model()
            self._construct_normalized_description()
            self._construct_equivalent_descriptions()

        # self.print_models_summary()

    def output(self):
        output = {
            "input_model": self.input_model,
        }
        if self.corrected_model:
            output["corrected_model"] = self.corrected_model
            output["corrected_description"] = model_to_string(self.corrected_model)
        output["normalized_description"] = self.normalized_description
        output["normalized_model"] = self.de_hgvs_model
        output["input_description"] = self.input_description

        if self.equivalent_descriptions is not None:
            output["equivalent_descriptions"] = self.equivalent_descriptions
        if self.errors:
            output["errors"] = self.errors
        if self.infos:
            output["infos"] = self.infos
        return output

    def print_models_summary(self):
        print("------")
        if self.input_description:
            print(self.input_description)

        if self.corrected_model:
            print("- Corrected model")
            print(model_to_string(self.corrected_model))
        else:
            print("- No corrected model")

        if self.internal_coordinates_model:
            print("- Internal coordinates model")
            print(model_to_string(self.internal_coordinates_model))
        else:
            print("- No internal coordinates model")

        if self.internal_indexing_model:
            print("- Internal indexing model")
            print(model_to_string(self.internal_indexing_model))
        else:
            print("- No internal_indexing_model")

        if self.delins_model:
            print("- Delins model")
            print(model_to_string(self.delins_model))
        else:
            print("- No delins model")

        if self.de_model:
            print("- De model")
            print(model_to_string(self.de_model))
        else:
            print("- No de_model")

        if self.de_hgvs_internal_indexing_model:
            print("- De hgvs internal indexing model")
            print(model_to_string(self.de_hgvs_internal_indexing_model))
        else:
            print("- De hgvs internal indexing model")

        if self.de_hgvs_coordinate_model:
            print("- De hgvs coordinate model")
            print(model_to_string(self.de_hgvs_coordinate_model))
        else:
            print("- No de_hgvs_coordinate_model")

        if self.de_hgvs_model:
            print("- De hgvs model")
            print(model_to_string(self.de_hgvs_model))
        else:
            print("- No de hgvs model")
        print("------")

    def get_reference_summary(self):
        return {
            "sequence_length": get_sequence_length(self.references, "reference"),
            "selector_ids": len(
                get_selectors_ids(self.references["reference"]["annotations"])
            ),
        }


def normalize(description_to_normalize):
    description = Description(description_to_normalize)
    description.normalize()
    return description.output()
