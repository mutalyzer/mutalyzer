import copy
import itertools

from Bio.Seq import Seq
from Bio.SeqUtils import seq1, seq3
from extractor import describe_dna
from mutalyzer_backtranslate import BackTranslate
from mutalyzer_hgvs_parser import to_model
from mutalyzer_hgvs_parser.exceptions import UnexpectedCharacter, UnexpectedEnd
from mutalyzer_mutator import mutate
from mutalyzer_mutator.util import reverse_complement

import mutalyzer.errors as errors
import mutalyzer.infos as infos

from .checker import (
    are_sorted,
    contains_insert_length,
    contains_uncertain_locations,
    is_overlap,
    splice_sites,
)
from .converter.extras import (
    convert_amino_acids,
    convert_reference_model,
    convert_selector_model,
)
from .converter.to_delins import to_delins, variants_to_delins
from .converter.to_hgvs_coordinates import (
    crossmap_to_hgvs_setup,
    point_to_hgvs,
    to_hgvs_locations,
)
from .converter.to_internal_coordinates import (
    crossmap_to_internal_setup,
    get_coordinate_system,
    point_to_internal,
    to_internal_coordinates,
)
from .converter.to_internal_indexing import to_internal_indexing
from .converter.to_rna import to_rna_reference_model, to_rna_sequences, to_rna_variants
from .converter.variants_de_to_hgvs import de_to_hgvs
from .description_model import (
    get_locations_min_max,
    get_reference_id,
    get_selector_id,
    model_to_string,
    yield_reference_ids,
    yield_reference_selector_ids,
    yield_reference_selector_ids_coordinate_system,
    yield_sub_model,
    yield_values,
)
from .protein import (
    get_protein_description,
    get_protein_sequence,
    in_frame_description,
)
from .reference import (
    ASSEMBLIES,
    ASSEMBLY_ALIASES,
    get_chromosome_accession,
    get_coordinate_system_from_reference,
    get_coordinate_system_from_selector_id,
    get_gene_selectors,
    get_gene_selectors_hgnc,
    get_internal_selector_model,
    get_only_selector_id,
    get_protein_selector_model,
    get_reference_id_from_model,
    get_reference_mol_type,
    get_selectors_ids,
    get_sequence_length,
    is_only_one_selector,
    is_selector_in_reference,
    overlap_min_max,
    retrieve_reference,
    slice_to_selector,
    yield_overlap_ids,
)
from .util import (
    check_errors,
    construct_sequence,
    get_end,
    get_location_length,
    get_start,
    get_submodel_by_path,
    is_dna,
    is_rna,
    point_in_insertion,
    reverse_path,
    set_by_path,
    slice_sequence,
    sort_variants,
)


class Description(object):
    def __init__(
        self,
        description=None,
        description_model=None,
        only_variants=False,
        sequence=None,
        stop_on_error=False,
    ):

        self.stop_on_errors = stop_on_error

        self.errors = []
        self.infos = []
        self.input_description = description if description else None
        self.input_model = description_model if description_model else {}
        self.only_variants = only_variants
        self.sequence = sequence
        self._check_input()

        self.corrected_model = copy.deepcopy(self.input_model)
        self.internal_coordinates_model = {}
        self.internal_indexing_model = {}
        self.delins_model = {}
        self.de_model = {}
        self.de_hgvs_internal_indexing_model = {}
        self.de_hgvs_coordinate_model = {}
        self.de_hgvs_model = {}
        self.normalized_description = None
        self.chromosomal_descriptions = None
        self.protein = None
        self.rna = None

        self.references = {}

        self.equivalent = []
        self.back_translated_descriptions = []

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
                errors.mismatch(model_description, self.input_description)

    def _add_error(self, error):
        self.errors.append(error)

    def _add_info(self, info):
        self.infos.append(info)

    def get_selector_id(self):
        if self.corrected_model:
            return get_selector_id(self.corrected_model)

    def get_selector_model(self):
        selector_id = self.get_selector_id()
        if self.references and selector_id:
            return get_internal_selector_model(
                self.references["reference"]["annotations"], selector_id, True
            )

    def is_inverted(self):
        selector_model = self.get_selector_model()
        if (
            selector_model
            and selector_model.get("location")
            and selector_model["location"].get("strand") == -1
        ):
            return True
        return False

    @check_errors
    def _convert_description_to_model(self):
        start_rule = "variants" if self.only_variants else "description"
        try:
            model = to_model(self.input_description, start_rule)
        except UnexpectedCharacter as e:
            self._add_error(errors.syntax_uc(e))
        except UnexpectedEnd as e:
            self._add_error(errors.syntax_ueof(e))
        else:
            if start_rule == "variants":
                self.input_model = {"variants": model}
            else:
                self.input_model = model

    def _set_main_reference(self):
        reference_id = get_reference_id(self.corrected_model)
        if reference_id in self.references:
            self.references["reference"] = self.references[reference_id]

    def _correct_reference_id(self, path, original_id, corrected_id):
        set_by_path(self.corrected_model, path, corrected_id)
        self._add_info(infos.corrected_reference_id(original_id, corrected_id, path))

    def _correct_lrg_reference_id(self, original_id, model, path):
        set_by_path(self.corrected_model, path[:-1], model)
        self._add_info(infos.corrected_lrg_reference(original_id, model, path))

    @staticmethod
    def _check_if_lrg_reference(reference_id):
        if reference_id.startswith("LRG_") and "t" in reference_id:
            lrg_split = reference_id.split("t")
            if len(lrg_split) == 2 and lrg_split[1].isdigit():
                return {
                    "id": lrg_split[0],
                    "selector": {"id": "t" + lrg_split[1]},
                }

    @check_errors
    def retrieve_references(self):
        """
        Populate the references
        """
        if not self.corrected_model:
            return
        if self.only_variants and self.sequence:
            self.references["reference"] = {"sequence": {"seq": self.sequence}}
        for reference_id, path in yield_reference_ids(self.corrected_model):
            reference_model = retrieve_reference(reference_id)[0]
            if reference_model is None:
                lrg = self._check_if_lrg_reference(reference_id)
                if lrg:
                    reference_model = retrieve_reference(lrg["id"])[0]
                    if reference_model:
                        self._correct_lrg_reference_id(reference_id, lrg, path)
                        reference_id = lrg["id"]
            if reference_model is None:
                self._add_error(errors.reference_not_retrieved(reference_id, [path]))
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
            infos.corrected_selector_id(
                original_id, corrected_id, correction_source, path
            )
        )

    def handle_selector_not_found(self, reference_id, selector_id, path):
        """
        Checks if the `selector_id` is either a gene name, HGNC gene id, or
        uses the legacy format (e.g., `SDHD_v1`) and updates the description
        model correspondingly, based on the provided tree path.

        :param reference_id: The reference id.
        :param selector_id: The selector id.
        :param path: Path in the description model tree.
        """
        gene_selectors = get_gene_selectors(selector_id, self.references[reference_id])
        if len(gene_selectors) == 1:
            self._correct_selector_id(path, selector_id, gene_selectors[0], "gene name")
            return
        elif len(gene_selectors) > 1:
            self._add_error(
                errors.selector_options(selector_id, "gene", gene_selectors, path)
            )
            return
        if "_v" in selector_id:
            gene_name = selector_id.split("_v")[0]
            gene_selectors = get_gene_selectors(
                gene_name, self.references[reference_id]
            )
            if len(gene_selectors) == 1:
                self._correct_selector_id(
                    path, selector_id, gene_selectors[0], "gene name"
                )
                return
            elif len(gene_selectors) > 1:
                self._add_error(
                    errors.selector_options(gene_name, "gene", gene_selectors, path)
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
                errors.selector_options(selector_id, "gene HGNC", gene_selectors, path)
            )
            return
        self._add_error(errors.no_selector_found(reference_id, selector_id, path))

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
            infos.corrected_coordinate_system(
                coordinate_system, correction_source, path
            )
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
        self._add_error(errors.no_coordinate_system(c_s_path))

    def _correct_selector_id_from_coordinate_system(self, r_id, r_id_path, selector_id):
        path = tuple(list(r_id_path[:-1]) + ["selector"])
        set_by_path(self.corrected_model, path, {"id": selector_id})
        if r_id != selector_id:
            self._add_info(
                infos.corrected_selector_id("", selector_id, "coordinate system", path)
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
                c_s_s = get_coordinate_system_from_selector_id(
                    self.references[r_id], s_id
                )
                if c_s_s != c_s and not (
                    (c_s_s == "c" and c_s in ["r", "p"])
                    or (c_s_s == "n" and c_s == "r")
                ):
                    self._add_error(
                        errors.coordinate_system_mismatch(c_s, s_id, c_s_s, c_s_path)
                    )
            else:
                r_c_s = get_coordinate_system_from_reference(self.references[r_id])
                if not ((r_c_s == c_s) and (c_s in ["g", "m", "p"])):
                    if is_only_one_selector(self.references[r_id]):
                        self._correct_selector_id_from_coordinate_system(
                            r_id,
                            r_path,
                            get_only_selector_id(self.references[r_id]),
                        )
                    elif r_c_s == "m" and c_s == "g":
                        self._correct_coordinate_system(
                            r_c_s, c_s_path, r_id + " reference"
                        )
                    else:
                        self._add_error(
                            errors.coordinate_system_mismatch(
                                c_s, r_id, r_c_s, c_s_path
                            )
                        )

    @check_errors
    def _check_selector_models(self):
        for reference_id, selector_id, path in yield_reference_selector_ids(
            self.corrected_model
        ):
            selector_model = get_internal_selector_model(
                self.references[reference_id]["annotations"], selector_id, True
            )
            if (
                self.corrected_model["coordinate_system"] in ["c", "p"]
                and selector_model.get("cds")
                and len(selector_model["cds"]) > 1
                and selector_model.get("exception")
                and selector_model.get("exception") == "ribosomal slippage"
            ):
                self._add_error(
                    errors.cds_slices(
                        f"Expcetion: \"{selector_model.get('exception')}\"."
                    )
                )
            if selector_model.get("whole_exon_transcript"):
                self._add_info(
                    infos.whole_transcript_exon(reference_id, selector_id, path)
                )

    @check_errors
    def _correct_variants_type(self):
        for i, v in enumerate(self.internal_indexing_model["variants"]):
            if v.get("type") == "substitution":
                if len(construct_sequence(v["inserted"], self.get_sequences())) > 1:
                    path = ["variants", i, "type"]
                    if self.is_inverted():
                        path = reverse_path(self.corrected_model, path)
                    set_by_path(self.corrected_model, path, "deletion_insertion")
                    set_by_path(
                        self.internal_coordinates_model, path, "deletion_insertion"
                    )
                    set_by_path(
                        self.internal_indexing_model, path, "deletion_insertion"
                    )
                    self._add_info(
                        infos.corrected_variant_type(
                            "substitution", "deletion insertion"
                        )
                    )

    @check_errors
    def _correct_chromosome_points(self):
        for point, path in yield_sub_model(
            self.corrected_model, ["location", "start", "end"], ["point"]
        ):
            if point["position"] in ["pter", "qter"]:
                ref_id = self.corrected_model["reference"]["id"]
                sel_id = self.get_selector_id()
                c_s = get_coordinate_system(self.corrected_model, self.references)
                for ins_or_del in ["inserted", "deleted"]:
                    if ins_or_del in path:
                        ins_or_del_model = get_submodel_by_path(
                            self.corrected_model,
                            path[: path.index(ins_or_del) + 2],
                        )
                        ins_or_del_ref_id = get_reference_id(ins_or_del_model)
                        ins_or_del_sel_id = get_selector_id(ins_or_del_model)
                        ins_or_del_c_s = ins_or_del_model.get("coordinate_system")
                        if ins_or_del_ref_id:
                            ref_id = ins_or_del_ref_id
                        if ins_or_del_sel_id:
                            sel_id = ins_or_del_sel_id
                        if ins_or_del_c_s:
                            c_s = ins_or_del_c_s
                sel_model = get_internal_selector_model(
                    self.references[ref_id]["annotations"], sel_id, True
                )
                crossmap_from = crossmap_to_hgvs_setup(c_s, sel_model, degenerate=True)
                if point["position"] == "pter":
                    to_correct = 0
                else:
                    to_correct = get_sequence_length(self.references, ref_id) - 1
                corrected = point_to_hgvs(
                    {"type": "point", "position": to_correct}, **crossmap_from
                )
                set_by_path(self.corrected_model, path, corrected)
                self._add_info(infos.corrected_point(point, corrected, path))

    @check_errors
    def _correct_points(self):
        c_s = get_coordinate_system(self.corrected_model, self.references)
        crossmap_to = crossmap_to_internal_setup(c_s, self.get_selector_model())
        crossmap_from = crossmap_to_hgvs_setup(
            c_s,
            self.get_selector_model(),
            degenerate=True,
        )

        for point, path in yield_sub_model(
            self.corrected_model, ["location", "start", "end"], ["point"]
        ):
            internal = point_to_internal(point, crossmap_to)
            corrected = point_to_hgvs(internal, **crossmap_from)
            if corrected != point:
                set_by_path(self.corrected_model, path, corrected)
                self._add_info(infos.corrected_point(point, corrected, path))

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
        if not are_sorted(self.delins_model["variants"]):
            self.delins_model["variants"] = sort_variants(self.delins_model["variants"])
            self._add_info(infos.sorted_variants())

    def get_sequences(self):
        """
        Retrieves a dictionary from the references with reference ids as
        keys and their corresponding sequences as values.
        """
        sequences = {k: self.references[k]["sequence"]["seq"] for k in self.references}
        return sequences

    @check_errors
    def _mutate(self):
        observed_sequence = mutate(self.get_sequences(), self.delins_model["variants"])
        self.references["observed"] = {"sequence": {"seq": observed_sequence}}

    @check_errors
    def _extract(self):
        self.de_model = {
            "variants": describe_dna(
                self.references["reference"]["sequence"]["seq"],
                self.references["observed"]["sequence"]["seq"],
            ),
        }
        if not self.only_variants:
            self.de_model.update(
                {
                    "reference": copy.deepcopy(
                        self.internal_indexing_model["reference"]
                    ),
                    "coordinate_system": "i",
                }
            )

    @check_errors
    def construct_de_hgvs_internal_indexing_model(self):
        if self.de_model:
            self.de_hgvs_internal_indexing_model = {
                "variants": de_to_hgvs(
                    self.de_model["variants"],
                    self.get_sequences(),
                ),
            }
            if not self.only_variants:
                self.de_hgvs_internal_indexing_model.update(
                    {
                        "reference": copy.deepcopy(
                            self.internal_indexing_model["reference"]
                        ),
                        "coordinate_system": "i",
                    }
                )
            if self.corrected_model.get("predicted"):
                self.de_hgvs_internal_indexing_model["predicted"] = True

    @check_errors
    def construct_de_hgvs_coordinates_model(self):
        if self.de_hgvs_internal_indexing_model:
            self.de_hgvs_model = to_hgvs_locations(
                self.de_hgvs_internal_indexing_model,
                self.references,
                self.corrected_model.get("coordinate_system"),
                get_selector_id(self.corrected_model),
                True,
            )

    def construct_normalized_description(self):
        if self.de_hgvs_model:
            if self.de_hgvs_model.get("coordinate_system") == "r":
                to_rna_sequences(self.de_hgvs_model)
            self.normalized_description = model_to_string(self.de_hgvs_model)

    def construct_equivalent(self, other=None, as_description=True):
        if self.only_variants:
            return
        if other is not None:
            from_model = other
        elif self.de_hgvs_internal_indexing_model:
            from_model = self.de_hgvs_internal_indexing_model
        else:
            return
        equivalent = {}

        if (
            get_coordinate_system_from_reference(self.references["reference"])
            == "g"
            != self.corrected_model["coordinate_system"]
            and self.corrected_model["coordinate_system"] != "r"
        ):
            converted_model = to_hgvs_locations(
                model=from_model,
                references=self.references,
                to_coordinate_system="g",
                to_selector_id=None,
                degenerate=True,
            )
            if as_description:
                equivalent["g"] = [model_to_string(converted_model)]
            else:
                equivalent["g"] = [converted_model]
            if equivalent:
                self.equivalent = equivalent

        l_min, l_max = get_locations_min_max(from_model)
        if not (l_min and l_max):
            return

        l_min, l_max = overlap_min_max(self.references["reference"], l_min, l_max)
        for selector in yield_overlap_ids(self.references["reference"], l_min, l_max):
            if selector["id"] != self.get_selector_id():
                converted_model = to_hgvs_locations(
                    model=from_model,
                    references=self.references,
                    to_coordinate_system=None,
                    to_selector_id=selector["id"],
                    degenerate=True,
                )
                c_s = converted_model["coordinate_system"]
                if not equivalent.get(c_s):
                    equivalent[c_s] = []

                if converted_model["coordinate_system"] == "c":
                    protein_selector_model = get_protein_selector_model(
                        self.references["reference"]["annotations"], selector["id"]
                    )
                    if protein_selector_model and as_description:
                        equivalent[c_s].append(
                            (
                                model_to_string(converted_model),
                                get_protein_description(
                                    variants_to_delins(from_model["variants"]),
                                    self.references,
                                    protein_selector_model,
                                )[0],
                            )
                        )
                    else:
                        if as_description:
                            equivalent[c_s].append(model_to_string(converted_model))
                        else:
                            equivalent[c_s].append(converted_model)
                else:
                    if as_description:
                        equivalent[c_s].append(model_to_string(converted_model))
                    else:
                        equivalent[c_s].append(converted_model)

        if equivalent:
            self.equivalent = equivalent

    @check_errors
    def construct_protein_description(self):
        if self.de_hgvs_model.get("coordinate_system") == "c" or (
            self.de_hgvs_model.get("coordinate_system") == "r"
            and get_coordinate_system_from_selector_id(
                self.references["reference"], self.get_selector_id()
            )
            == "c"
        ):
            if self.rna and self.rna.get("errors"):
                self.protein = {"errors": self.rna["errors"]}
            else:
                protein_selector_model = get_protein_selector_model(
                    self.references["reference"]["annotations"],
                    get_selector_id(self.de_hgvs_model),
                )
                if protein_selector_model:
                    self.protein = dict(
                        zip(
                            ["description", "reference", "predicted"],
                            get_protein_description(
                                variants_to_delins(
                                    self.de_hgvs_internal_indexing_model["variants"]
                                ),
                                self.references,
                                protein_selector_model,
                            ),
                        )
                    )

    @check_errors
    def construct_rna_description(self):
        if self.de_hgvs_model.get("coordinate_system") in ["c", "n"]:
            self.rna = {}
            errors_splice, infos_splice = splice_sites(
                variants_to_delins(self.de_hgvs_internal_indexing_model["variants"]),
                self.get_sequences(),
                self.get_selector_model(),
            )
            if infos_splice:
                self.rna["infos"] = infos_splice
            if errors_splice:
                self.rna["errors"] = errors_splice
                return

            rna_variants_coordinate = to_rna_variants(
                variants_to_delins(self.de_hgvs_internal_indexing_model["variants"]),
                self.get_sequences(),
                self.get_selector_model(),
            )
            rna_reference_model = to_rna_reference_model(
                self.references["reference"], self.get_selector_id()
            )
            rna_references = {
                get_reference_id(self.corrected_model): rna_reference_model,
                "reference": rna_reference_model,
            }
            rna_variants_coordinate = de_to_hgvs(
                rna_variants_coordinate,
                {
                    k: str(
                        Seq(rna_references[k]["sequence"]["seq"]).transcribe().lower()
                    )
                    for k in rna_references
                },
            )
            to_rna_sequences(rna_variants_coordinate)
            rna_model = to_hgvs_locations(
                {
                    "reference": self.de_hgvs_internal_indexing_model["reference"],
                    "coordinate_system": "i",
                    "variants": rna_variants_coordinate,
                },
                rna_references,
                self.corrected_model["coordinate_system"],
                get_selector_id(self.corrected_model),
                True,
            )
            rna_model["coordinate_system"] = "r"
            rna_model["predicted"] = True
            self.rna = {"description": model_to_string(rna_model)}

    def _get_reference_id(self, model, path):
        ref_id = get_reference_id(model)
        for ins_or_del in ["inserted", "deleted"]:
            if ins_or_del in path:
                ins_or_del_ref = get_reference_id(
                    get_submodel_by_path(
                        model,
                        path[: path.index(ins_or_del) + 2],
                    )
                )
                if ins_or_del_ref:
                    ref_id = ins_or_del_ref
        return ref_id

    def _check_location_boundaries(self):
        for point, path in yield_sub_model(
            self.internal_coordinates_model, ["location", "start", "end"], ["point"]
        ):
            ref_id = "reference"
            for ins_or_del in ["inserted", "deleted"]:
                if ins_or_del in path:
                    ins_or_del_ref = get_reference_id(
                        get_submodel_by_path(
                            self.internal_coordinates_model,
                            path[: path.index(ins_or_del) + 2],
                        )
                    )
                    if ins_or_del_ref:
                        ref_id = ins_or_del_ref
            len_seq = get_sequence_length(self.references, ref_id)
            if point_in_insertion(self.internal_coordinates_model, path):
                left_boundary = -1
                right_boundary = len_seq + 1
            else:
                left_boundary = 0
                right_boundary = len_seq
            if not point.get("uncertain") and point.get("position"):
                point = point["position"]
                if right_boundary <= point:
                    if self.is_inverted():
                        path = reverse_path(self.internal_coordinates_model, path)
                    self._add_error(
                        errors.out_of_boundary_greater(
                            get_submodel_by_path(self.corrected_model, path),
                            point - right_boundary + 1,
                            right_boundary,
                            path,
                        )
                    )
                elif point < left_boundary:
                    if self.is_inverted():
                        path = reverse_path(self.internal_coordinates_model, path)
                    self._add_error(
                        errors.out_of_boundary_lesser(
                            get_submodel_by_path(self.corrected_model, path),
                            -point,
                            path,
                        )
                    )

    def _check_location_range(self):
        for location, path in yield_sub_model(
            self.internal_coordinates_model, ["location"], ["range"]
        ):
            if (
                not location["start"].get("uncertain")
                and not location["end"].get("uncertain")
            ) and get_start(location) > get_end(location):
                self._add_error(errors.range_reversed(location, path))

    def _check_genomic_point(self, point, path):
        if point.get("offset") or point.get("outside_cds"):
            c_s = self.corrected_model.get("coordinate_system")
            for ins_or_del in ["inserted", "deleted"]:
                if ins_or_del in path:
                    ins_or_del_c_s = get_submodel_by_path(
                        self.corrected_model,
                        path[: path.index(ins_or_del) + 2],
                    ).get("coordinate_system")
                    if ins_or_del_c_s:
                        c_s = ins_or_del_c_s
            if c_s == "g":
                if point.get("offset"):
                    self._add_error(errors.offset(point, path))
                if point.get("outside_cds"):
                    self._add_error(errors.outside_cds(point, path))

    def _check_intronic_point(self, point, path):
        if point.get("offset"):
            ref_id = self.corrected_model["reference"]["id"]
            for ins_or_del in ["inserted", "deleted"]:
                if ins_or_del in path:
                    ins_or_del_ref_id = get_reference_id(
                        get_submodel_by_path(
                            self.corrected_model,
                            path[: path.index(ins_or_del) + 2],
                        )
                    )
                    if ins_or_del_ref_id:
                        ref_id = ins_or_del_ref_id
            if get_reference_mol_type(self.references[ref_id]) in [
                "mRNA",
                "ncRNA",
                "transcribed RNA",
            ]:
                if point.get("offset"):
                    # TODO: find the actual NM(NC) description
                    self._add_error(errors.intronic(point, path))

    def _check_location_extras(self):
        for point, path in yield_sub_model(
            self.corrected_model, ["location", "start", "end"], ["point"]
        ):
            self._check_genomic_point(point, path)
            self._check_intronic_point(point, path)

    def _check_insertion_location(self, path):
        """
        Checks if the insertion location range is according to HGVS.
        Incorrect descriptions example:
            NM_000143.3:c.40_50insATC
        Notes:
            - One could argue to correct it to a delins. This can be done as
            a proposal in the web interface.
        :param path: Model path towards the variant ["variants", #].
        """
        v = self.internal_coordinates_model["variants"][path[1]]
        v_r = self.input_model["variants"][path[1]]

        if v["location"]["type"] == "point" and not v["location"].get("uncertain"):
            self._add_error(errors.insertion_range(v_r["location"], path))

        if v["location"]["type"] == "range" and not v["location"].get("uncertain"):
            if abs(get_start(v) - get_end(v)) != 1:
                self._add_error(errors.insertion_range(v_r["location"], path))

    def _check_repeat_dbsnp(self, path, v, v_i):
        ref_seq = self.references["reference"]["sequence"]["seq"]
        inserted = v_i["inserted"][0]
        if self.is_inverted():
            repeat_unit = reverse_complement(inserted["sequence"])
            point = get_start(v_i)
            while ref_seq[point - len(repeat_unit) + 1 : point + 1] == repeat_unit:
                point -= len(repeat_unit)
            if point + 1 < get_start(v_i):
                v_i["location"]["start"]["position"] = point + 1
            else:
                self._add_error(
                    errors.repeat_sequences_mismatch("", inserted["sequence"], path)
                )
        else:
            inserted = v_i["inserted"][0]
            repeat_unit = inserted["sequence"]
            point = get_start(v_i)
            while ref_seq[point : point + len(repeat_unit)] == repeat_unit:
                point += len(repeat_unit)
            if point > get_start(v_i):
                v_i["location"]["end"]["position"] = point
            else:
                self._add_error(errors.repeat_sequences_mismatch("", repeat_unit, path))

    def _check_repeat(self, path):
        if self.is_inverted():
            path_i = reverse_path(self.internal_indexing_model, path)
        else:
            path_i = path
        v = self.input_model["variants"][path_i[1]]
        v_i = self.internal_indexing_model["variants"][path[1]]
        if v_i.get("inserted") and len(v_i.get("inserted")) == 1:
            inserted = v_i["inserted"][0]
            if v["location"]["type"] == "point" and len(inserted) > 1:
                self._check_repeat_dbsnp(path, v, v_i)
                return
            if inserted.get("sequence") and inserted.get("source") == "description":
                repeat_seq = inserted["sequence"]
            # TODO: get the sequence from a reference slice
            else:
                self._add_error(errors.repeat_not_supported(v, path))
                return
            ref_seq = self.references["reference"]["sequence"]["seq"][
                get_start(v_i) : get_end(v_i)
            ]
            if self.is_inverted():
                ref_seq = reverse_complement(ref_seq)

            if len(ref_seq) % len(repeat_seq) != 0:
                self._add_error(errors.repeat_reference_sequence_length(path))
            elif (len(ref_seq) // len(repeat_seq)) * repeat_seq != ref_seq:
                self._add_error(
                    errors.repeat_sequences_mismatch(ref_seq, repeat_seq, path)
                )
        else:
            # TODO: Convert to delins and switch to warning?
            self._add_error(errors.repeat_not_supported(v, path))

    def _check_superfluous(self, path):
        """
        Checks if any additional provided sequence/length is correct.
        Incorrect descriptions examples:
            NM_000143.3:c.45delT
            NG_012337.1:g.274ATT>T
            NM_000143.3:c.45dupT
            NM_000143.3:c.45dup4
        :param path: Model path towards "deleted" or "inserted".
        """
        v_i = self.internal_indexing_model["variants"][path[1]]
        ins_or_del = path[-1]
        sequences = self.get_sequences()
        if (
            len(v_i[ins_or_del]) == 1
            and v_i[ins_or_del][0].get("length")
            and v_i[ins_or_del][0]["length"].get("value")
        ):
            len_del = v_i[ins_or_del][0]["length"].get("value")
            len_loc = get_location_length(v_i["location"])
            if len_loc != len_del:
                self._add_error(errors.length_mismatch(len_loc, len_del, path))
        elif get_start(v_i["location"]) >= 0 and get_end(v_i["location"]) < len(
            sequences["reference"]
        ):
            seq_ref = slice_sequence(v_i["location"], sequences["reference"])
            seq_del = construct_sequence(v_i[ins_or_del], sequences)
            if self.corrected_model.get("coordinate_system") == "r":
                seq_del = str(Seq(seq_del).transcribe().lower())
                seq_ref = str(Seq(seq_ref).transcribe().lower())
            if seq_del and seq_ref and seq_del != seq_ref:
                if self.is_inverted():
                    seq_del = reverse_complement(seq_del)
                    seq_ref = reverse_complement(seq_ref)
                    path = reverse_path(self.internal_coordinates_model, path)
                self._add_error(errors.sequence_mismatch(seq_ref, seq_del, path))

    @check_errors
    def _check_and_correct_sequences(self):
        for seq, path in yield_sub_model(
            self.corrected_model, ["sequence", "amino_acid"]
        ):
            if self.corrected_model.get("coordinate_system") in ["g", "c", "n", None]:
                if seq.upper() != seq:
                    self._add_info(infos.corrected_sequence(seq, seq.upper()))
                    set_by_path(self.corrected_model, path, seq.upper())
                    if self.is_inverted():
                        path = reverse_path(self.corrected_model, path)
                    set_by_path(self.internal_coordinates_model, path, seq.upper())
                    set_by_path(self.internal_indexing_model, path, seq.upper())
                if not is_dna(seq.upper()):
                    self._add_error(errors.no_dna(seq.upper(), path))
            elif self.corrected_model.get("coordinate_system") == "r":
                if is_dna(seq):
                    seq_rna = str(Seq(seq).transcribe()).lower()
                    self._add_info(infos.corrected_sequence(seq, seq_rna))
                    seq_dna = str(Seq(seq_rna).back_transcribe()).upper()
                    set_by_path(self.corrected_model, path, seq_dna)
                    if self.is_inverted():
                        path = reverse_path(self.corrected_model, path)
                    set_by_path(self.internal_coordinates_model, path, seq_dna)
                    set_by_path(self.internal_indexing_model, path, seq_dna)
                else:
                    seq_lower = seq.lower()
                    if seq_lower != seq:
                        self._add_info(infos.corrected_sequence(seq, seq_lower))
                    if not is_rna(seq_lower):
                        self._add_error(errors.no_rna(seq_lower, path))
                    seq_lower = str(Seq(seq).back_transcribe()).upper()
                    # set_by_path(self.corrected_model, path, seq_lower)
                    if self.is_inverted():
                        path = reverse_path(self.corrected_model, path)
                    set_by_path(self.internal_coordinates_model, path, seq_lower)
                    set_by_path(self.internal_indexing_model, path, seq_lower)

    @check_errors
    def _check_location_amino_acids(self):
        sequences = self.get_sequences()
        for point, path in yield_sub_model(
            self.internal_coordinates_model, ["location", "start", "end"], ["point"]
        ):
            ref_id = self._get_reference_id(self.internal_indexing_model, path)
            if (
                point.get("amino_acid")
                and point.get("position")
                and sequences[ref_id][point["position"]] != point["amino_acid"]
            ):
                self._add_error(
                    errors.amino_acid_mismatch(
                        point["amino_acid"], sequences[ref_id][point["position"]], path
                    )
                )

    def _insertions_same_location(self):
        insertions = {}
        variants = self.internal_coordinates_model["variants"]
        for i, v in enumerate(variants):
            if v.get("type") == "insertion":
                s_e = (get_start(v), get_end(v))
                if s_e not in insertions:
                    insertions[s_e] = []
                insertions[s_e].append(i)
        for insertion in insertions.values():
            if len(insertion) > 1:
                return True
        return False

    @check_errors
    def check(self):
        self._check_location_boundaries()
        self._check_location_range()
        self._check_location_amino_acids()
        for i, v in enumerate(self.internal_coordinates_model["variants"]):
            if v.get("deleted"):
                self._check_superfluous(["variants", i, "deleted"])

            if v.get("type") in ["duplication", "inversion"] and v.get("inserted"):
                self._check_superfluous(["variants", i, "inserted"])

            if v.get("type") == "insertion":
                self._check_insertion_location(["variants", i])

            if v.get("type") == "repeat":
                self._check_repeat(["variants", i])
        if (
            is_overlap(self.internal_indexing_model["variants"])
            or self._insertions_same_location()
        ):
            self._add_error(errors.overlap())

    @check_errors
    def pre_conversion_checks(self):
        self._check_selectors_in_references()
        self._check_coordinate_systems()
        self._check_coordinate_system_consistency()
        self._check_selector_models()
        self._check_location_extras()
        if contains_uncertain_locations(self.corrected_model):
            self._add_error(errors.uncertain())
        if contains_insert_length(self.corrected_model):
            self._add_error(errors.inserted_length())

    def assembly_checks(self):
        r_id = get_reference_id(self.input_model)
        if r_id is not None:
            assembly = r_id.upper()
        else:
            return
        if assembly in ASSEMBLY_ALIASES:
            assembly = ASSEMBLY_ALIASES[assembly]
        elif assembly not in ASSEMBLIES:
            return

        s_id = get_selector_id(self.input_model)
        if s_id is not None:
            if s_id.upper().startswith("CHR"):
                chr_number = s_id.upper().split("CHR")[1]
            else:
                chr_number = s_id.upper()
        else:
            return

        if chr_number in ASSEMBLIES[assembly]:
            chromosome_id = ASSEMBLIES[assembly][chr_number]
            self.corrected_model["reference"] = self.corrected_model["reference"][
                "selector"
            ]
            self.corrected_model["reference"]["id"] = chromosome_id
            self._add_info(infos.assembly_chromosome_to_id(r_id, s_id, chromosome_id))

    def to_internal_indexing_model(self):
        self._construct_internal_coordinate_model()
        self._construct_internal_indexing_model()

    def remove_superfluous_selector(self):
        if (
            not self.only_variants
            and self.de_hgvs_model
            and self.de_hgvs_model["reference"].get("selector")
            and self.de_hgvs_model["reference"]["selector"]["id"]
            == self.de_hgvs_model["reference"]["id"]
        ):
            self.de_hgvs_model["reference"].pop("selector")

    @check_errors
    def only_equals(self):
        for variant in self.internal_coordinates_model["variants"]:
            if variant.get("type") != "equal":
                return False
        return True

    @check_errors
    def no_operation(self):
        if self.internal_coordinates_model.get("variants") is None:
            return True
        for variant in self.internal_coordinates_model["variants"]:
            if variant.get("type") is not None:
                return False
        return True

    @check_errors
    def _rna(self):
        if self.corrected_model.get("coordinate_system") == "r":
            errors_splice, infos_splice = splice_sites(
                self.internal_indexing_model["variants"],
                self.get_sequences(),
                self.get_selector_model(),
            )
            self.infos += infos_splice
            self.errors += errors_splice
            if errors_splice:
                return

            self.internal_indexing_model["variants"] = to_rna_variants(
                self.internal_indexing_model["variants"],
                self.get_sequences(),
                self.get_selector_model(),
            )
            rna_reference_model = to_rna_reference_model(
                self.references["reference"], self.get_selector_id()
            )
            # self.delins_model["variants"] = variants
            self.references = {
                get_reference_id(self.corrected_model): rna_reference_model,
                "reference": rna_reference_model,
            }

    def _check_amino_acids(self):
        for sequence, path in yield_values(
            self.corrected_model, ["sequence", "amino_acid"]
        ):
            seq_1a = str(seq1(sequence))
            seq_3a = str(seq3(sequence))
            if not ((sequence == str(seq3(seq_1a))) != (sequence == str(seq1(seq_3a)))):
                self._add_error(
                    {
                        "code": "EAA",
                        "details": f"Sequence '{sequence}' is a mix of 1 and 3 letter amino acids.",
                    }
                )

    def _check_supported_variants(self, model):
        for i, variant in enumerate(model["variants"]):
            if variant.get("type") in ["frame_shift", "extension"]:
                self._add_error(
                    errors.variant_not_supported(
                        variant, variant["type"], ["variants", i]
                    )
                )

    @check_errors
    def _convert_amino_acids(self):
        convert_amino_acids(self.internal_indexing_model, "1a")
        convert_amino_acids(self.internal_coordinates_model, "1a")

    def _back_translate(self):
        if not self.get_selector_id():
            return
        cds_seq = slice_to_selector(
            retrieve_reference(get_reference_id(self.corrected_model))[0],
            self.get_selector_id(),
            True,
            True,
        )
        bt = BackTranslate()
        translated_vars = []
        for variant in self.internal_indexing_model["variants"]:
            if variant.get("type") == "substitution":
                cds_start = get_start(variant) * 3
                cds_end = get_end(variant) * 3
                bt_options = bt.with_dna(
                    cds_seq[cds_start:cds_end],
                    variant["inserted"][0]["sequence"],
                )
                dna_var_options = []
                for offset in bt_options:
                    for v in bt_options[offset]:
                        dna_var_options.append(
                            "{}{}>{}".format(cds_start + offset + 1, v[0], v[1])
                        )
                translated_vars.append(dna_var_options)
            else:
                # TODO: Add error message.
                return []
        if self.get_selector_id() == get_reference_id(self.corrected_model):
            reference = "{}".format(get_reference_id(self.corrected_model))
        else:
            reference = "{}({})".format(
                get_reference_id(self.corrected_model),
                self.get_selector_model()["mrna_id"],
            )
        bt_descriptions = []
        for t in itertools.product(*translated_vars):
            if len(t) > 1:
                bt_descriptions.append("{}:c.([{}])".format(reference, ";".join(t)))
            elif len(t) == 1:
                bt_descriptions.append("{}:c.({})".format(reference, t[0]))
        self.back_translated_descriptions = bt_descriptions

    def normalize_protein(self):
        reference_id = self.corrected_model["reference"]["id"]
        self._check_amino_acids()
        if self.errors:
            return
        if self.get_selector_id():
            # Convert references to protein model
            p_seq = get_protein_sequence(
                self.references["reference"], self.get_selector_model()
            )
            self.references = copy.deepcopy(self.references)
            self.references["reference"]["sequence"]["seq"] = p_seq
            self.references[reference_id]["sequence"]["seq"] = p_seq

        self.to_internal_indexing_model()
        self._convert_amino_acids()
        self._correct_variants_type()
        self._check_and_correct_sequences()
        self.check()
        self._check_supported_variants(self.corrected_model)

        if self.errors:
            return

        if self.only_equals() or self.no_operation():
            self.de_hgvs_model = copy.deepcopy(self.corrected_model)
            convert_amino_acids(self.de_hgvs_model, "1a")

            convert_amino_acids(self.de_hgvs_model, "3a")
            self.references["observed"] = {
                "sequence": {"seq": self.references["reference"]["sequence"]["seq"]}
            }
        else:
            observed_sequence = mutate(
                self.get_sequences(),
                to_delins(self.internal_indexing_model)["variants"],
            )
            self.references["observed"] = {"sequence": {"seq": observed_sequence}}
            selector_id = (
                "({})".format(self.get_selector_id()) if self.get_selector_id() else ""
            )
            p_description = "{}{}:{}".format(
                reference_id,
                selector_id,
                in_frame_description(
                    self.references["reference"]["sequence"]["seq"], observed_sequence
                )[0],
            )
            self.de_hgvs_model = to_model(p_description)
        self.normalized_description = model_to_string(self.de_hgvs_model)
        equivalent_1a_model = copy.deepcopy(self.de_hgvs_model)
        convert_amino_acids(equivalent_1a_model, "1a")
        self.equivalent = {"p": [model_to_string(equivalent_1a_model)]}
        self._back_translate()

    def to_delins(self):
        self.retrieve_references()
        self.pre_conversion_checks()

        self._correct_chromosome_points()
        self.to_internal_indexing_model()
        self._correct_variants_type()
        self._correct_points()
        self._check_and_correct_sequences()

        self.check()
        self._rna()
        self._construct_delins_model()

    @check_errors
    def get_chromosomal_description(self):
        # TODO: Add tests.
        if not self.references or self.only_variants:
            return
        ref_id = get_reference_id(self.corrected_model)
        if (
            ref_id
            and get_reference_mol_type(self.references[ref_id]) == "mRNA"
            and self.corrected_model["coordinate_system"] == "c"
            and (ref_id.startswith("NM_") or ref_id.startswith("XM_"))
        ):
            self._add_info(infos.mrna_genomic_tip())
        elif (
            ref_id
            and get_reference_mol_type(self.references[ref_id]) == "mRNA"
            and self.corrected_model["coordinate_system"] == "r"
            and (ref_id.startswith("NM_") or ref_id.startswith("XM_"))
        ):
            return

        chromosome_accessions = get_chromosome_accession(
            ref_id, self.references["reference"]
        )
        if chromosome_accessions:
            chromosomal_descriptions = []
            for assembly, chromosome_accession in chromosome_accessions:
                chromosome_model = retrieve_reference(chromosome_accession)[0]
                if chromosome_model:
                    selector_ids = get_selectors_ids(
                        chromosome_model["annotations"], "c"
                    )
                    ref_id_accession = ref_id.split(".")[0]
                    matched_selector_ids = [
                        i for i in selector_ids if i.startswith(ref_id_accession)
                    ]
                    if not matched_selector_ids:
                        self._add_info(infos.no_selector(chromosome_accession, ref_id))
                        continue
                    if ref_id not in matched_selector_ids:
                        self._add_info(
                            infos.other_versions(
                                chromosome_accession, ref_id, matched_selector_ids
                            )
                        )
                        continue

                    selector_id = sorted(matched_selector_ids)[-1]

                    to_reference_model = convert_reference_model(
                        chromosome_model, selector_id, "transcript"
                    )

                    ref_seq_from = self.get_sequences()["reference"]
                    obs_seq = self.get_sequences()["observed"]
                    ref_seq_to = to_reference_model["sequence"]["seq"]

                    to_inverted = get_internal_selector_model(
                        chromosome_model["annotations"], selector_id, True
                    ).get("inverted")
                    if to_inverted:
                        obs_seq = reverse_complement(obs_seq)
                    variants = de_to_hgvs(
                        describe_dna(ref_seq_to, obs_seq),
                        {"reference": ref_seq_to, "observed": obs_seq},
                    )

                    if (not to_inverted and ref_seq_to != ref_seq_from) or (
                        to_inverted and reverse_complement(ref_seq_to) != ref_seq_from
                    ):
                        self._add_info(
                            infos.mrna_genomic_difference(ref_id, chromosome_accession)
                        )
                    else:
                        variants_model = to_hgvs_locations(
                            {
                                "reference": {
                                    "id": chromosome_accession,
                                    "selector": {"id": selector_id},
                                },
                                "coordinate_system": "i",
                                "variants": variants,
                            },
                            {
                                "reference": to_reference_model,
                                chromosome_model["annotations"][
                                    "id"
                                ]: to_reference_model,
                            },
                            get_coordinate_system_from_selector_id(
                                chromosome_model, selector_id
                            ),
                            selector_id,
                            True,
                        )
                        chromosomal_descriptions.append(
                            {
                                "assembly": assembly,
                                "description": model_to_string(variants_model),
                            }
                        )

            if chromosomal_descriptions:
                self.chromosomal_descriptions = chromosomal_descriptions

    def normalize_only_equals_or_no_operation(self):
        self.de_hgvs_internal_indexing_model = self.internal_indexing_model
        self.references["observed"] = {
            "sequence": {"seq": self.references["reference"]["sequence"]["seq"]}
        }
        self.construct_de_hgvs_coordinates_model()
        self.construct_normalized_description()
        self.construct_equivalent()

    def normalize(self):
        self.assembly_checks()
        self.retrieve_references()
        self.pre_conversion_checks()

        if self.corrected_model.get("type") == "description_protein":
            self.normalize_protein()
        else:
            self._correct_chromosome_points()
            self.to_internal_indexing_model()
            self._correct_variants_type()
            self._correct_points()
            self._check_and_correct_sequences()

            self.check()
            self._rna()
            self._construct_delins_model()
            if self.only_equals() or self.no_operation():
                self.normalize_only_equals_or_no_operation()
            else:
                self._mutate()
                self._extract()
                self.construct_de_hgvs_internal_indexing_model()
                self.construct_de_hgvs_coordinates_model()
                self.construct_normalized_description()
                self.construct_rna_description()
                self.construct_protein_description()
                self.construct_equivalent()
                self.get_chromosomal_description()
            self.remove_superfluous_selector()

        # self.print_models_summary()

    def output(self):
        output = {}
        if self.input_description:
            output["input_description"] = self.input_description
        if self.input_model:
            output["input_model"] = self.input_model
        if self.only_variants:
            output["only_variants"] = self.only_variants
            if self.sequence:
                output["sequence"] = self.sequence
        if self.corrected_model:
            output["corrected_model"] = self.corrected_model
            output["corrected_description"] = model_to_string(self.corrected_model)
        if self.normalized_description:
            output["normalized_description"] = self.normalized_description
        if self.de_hgvs_model:
            output["normalized_model"] = self.de_hgvs_model
        if self.chromosomal_descriptions:
            output["chromosomal_descriptions"] = self.chromosomal_descriptions
        if self.protein:
            output["protein"] = self.protein
        if self.rna:
            output["rna"] = self.rna
        if self.equivalent:
            output["equivalent_descriptions"] = self.equivalent
        if self.errors:
            output["errors"] = self.errors
        if self.infos:
            output["infos"] = self.infos
        if self.get_selector_model():
            output["selector_short"] = convert_selector_model(self.get_selector_model())
        if self.back_translated_descriptions:
            output["back_translated_descriptions"] = self.back_translated_descriptions
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
        if self.errors:
            print("---\nErrors:")
            print(self.errors)
        if self.infos:
            print("---\nInfos:")
            print(self.infos)
        print("------")

    def get_reference_summary(self):
        return {
            "sequence_length": get_sequence_length(self.references, "reference"),
            "selector_ids": len(
                get_selectors_ids(self.references["reference"]["annotations"])
            ),
        }

    def __str__(self):
        return self.input_description


def normalize(description_to_normalize):
    description = Description(description_to_normalize)
    description.normalize()
    return description.output()
