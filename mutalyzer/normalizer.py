"""Interfaces to obtain the normalized (canonical) variant representations or
the delins model of an input description."""

import itertools

from algebra import Variant
from algebra.extractor import extract as extract_variants
from algebra.extractor import local_supremal
from algebra.extractor import to_hgvs as to_hgvs_experimental
from algebra.lcs.lcs_graph import trim
from algebra.utils import to_dot
from mutalyzer_hgvs_parser import to_model

from mutalyzer.util import get_inserted_sequence

from .algebra import algebra_variant_to_delins, algebra_variant_to_name_model
from .converter.to_delins import to_delins, variants_to_delins
from .converter.to_hgvs_coordinates import (
    crossmap_to_hgvs_setup,
    initialize_hgvs_model,
    point_to_hgvs,
)
from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .description import Description
from .description_model import (
    get_selector_id,
    model_to_string,
    yield_point_locations_for_main_reference,
)
from . import infos
from .protein import get_protein_description, has_recoding_translation_exception
from .reference import get_protein_selector_model
from .rna import dna_to_rna, rna_to_dna
from .util import (
    construct_sequence,
    create_exact_point_model,
    create_exact_range_model,
    get_end,
    get_start,
    reverse_complement,
    set_by_path,
)
from .viewer import view_delins


def to_hgvs_dict(variants, ref_seq, forward_strand=True):
    """Algebra based experimental version of HGVS serialization with support for
    tandem repeats and complex variants."""
    def var_dict(var_type, start, end=None, inserted=None, repeat_number=None, del_seq=None):
        output = {
            "location": to_hgvs_position(start, end),
            "type": var_type,
            "source": "reference",
        }
        if isinstance(inserted, list):
            output["inserted"] = inserted
        else:
            if inserted:
                output["inserted"] = [{"sequence": inserted, "source": "description"}]
            if inserted and repeat_number:
                output["inserted"][0]["repeat_number"] = {"type": "point", "value": repeat_number}
        if del_seq:
            output["deleted"] = [{"sequence": del_seq, "source": "description"}]
        return output

    def repeats(word):
        length = 0
        idx = 1
        lps = [0] * len(word)
        while idx < len(word):
            if word[idx] == word[length]:
                length += 1
                lps[idx] = length
                idx += 1
            elif length != 0:
                length = lps[length - 1]
            else:
                lps[idx] = 0
                idx += 1

        pattern = len(word) - length
        if pattern == 0:
            return "", 0, 0
        return word[:pattern], len(word) // pattern, len(word) % pattern

    def to_hgvs_position(start, end=None):
        if end is None or end - start == 1:
            return create_exact_point_model(start + 1)
        if start == end:
            return create_exact_range_model(start, start + 1)
        return create_exact_range_model(start + 1, end)

    def other(variant):
        if variant.end - variant.start == 0:
            if not variant.sequence:
                return "="
            return var_dict("insertion", variant.start, variant.start, variant.sequence)
            # return f"{variant.start}_{variant.start + 1}ins{variant.sequence}"

        deleted = ""
        substitution = ref_seq[variant.start:variant.end]

        if variant.end - variant.start == 1:
            if not variant.sequence:
                return var_dict("deletion", variant.start)
                # return f"{variant.start + 1}del{deleted}"
            if len(variant.sequence) == 1:
                return var_dict("substitution", variant.start, del_seq=substitution, inserted=variant.sequence)
                # return f"{variant.start + 1}{substitution}>{variant.sequence}"
            return var_dict("deletion_insertion", variant.start, del_seq=deleted, inserted=variant.sequence)
            # return f"{variant.start + 1}del{deleted}ins{variant.sequence}"

        if not variant.sequence:
            return var_dict("deletion", variant.start, variant.end, del_seq=deleted)
            # return f"{variant.start + 1}_{variant.end}del{deleted}"

        return var_dict("deletion_insertion", variant.start, variant.end, del_seq=deleted, inserted=variant.sequence)
        # return f"{variant.start + 1}_{variant.end}del{deleted}ins{variant.sequence}"

    def hgvs(variant):
        inserted_unit, inserted_number, inserted_remainder = repeats(variant.sequence)
        deleted = ref_seq[variant.start:variant.end]
        deleted_unit, deleted_number, deleted_remainder = repeats(deleted)

        # Select a non-minimal repeat unit if reference and observed are
        # in agreement.
        diff = len(inserted_unit) - len(deleted_unit)
        if diff < 0 and deleted_unit == variant.sequence[:len(inserted_unit) - diff]:
            inserted_unit = deleted_unit
            inserted_number = 1
            inserted_remainder = deleted_remainder
        elif diff > 0 and inserted_unit == deleted[:len(deleted_unit) + diff]:
            deleted_unit = inserted_unit
            deleted_number = 1
            deleted_remainder = inserted_remainder

        # Repeat structure
        if deleted_unit == inserted_unit:
            if deleted_number == inserted_number:
                raise ValueError("empty variant")

            # Duplication
            if deleted_number == 1 and inserted_number == 2:
                return var_dict(
                    "duplication",
                    variant.start + inserted_remainder,
                    variant.start + inserted_remainder + len(inserted_unit),
                )
                # return f"{to_hgvs_position(variant.start + inserted_remainder, variant.start + inserted_remainder + len(inserted_unit))}dup"

            # shift 3'
            assert deleted_remainder == inserted_remainder
            inserted_unit = variant.sequence[inserted_remainder:inserted_remainder + len(inserted_unit)]

            return var_dict("repeat", variant.start + deleted_remainder, variant.end, inserted_unit, inserted_number)
            # return f"{to_hgvs_position(variant.start + deleted_remainder, variant.end)}{inserted_unit}[{inserted_number}]"

        # Prefix and suffix trimming
        start, end = trim(deleted, variant.sequence)
        trimmed = Variant(variant.start + start, variant.end - end, variant.sequence[start:len(variant.sequence) - end])

        # Inversion
        if len(trimmed.sequence) > 1 and trimmed.sequence == reverse_complement(ref_seq[trimmed.start:trimmed.end]):
            return var_dict("inversion", trimmed.start, trimmed.end)
            # return f"{to_hgvs_position(trimmed.start, trimmed.end)}inv"

        # Deletion/insertion with repeated insertion
        inserted_unit, inserted_number, inserted_remainder = repeats(trimmed.sequence)
        if inserted_number > 1:
            suffix = [{"sequence": inserted_unit, "source": "description", "repeat_number": {"type": "point", "value": inserted_number}}]
            # suffix = f"{inserted_unit}[{inserted_number}]"
            if inserted_remainder:
                suffix = [suffix[0], {"sequence": inserted_unit[:inserted_remainder], "source": "description"}]
                # suffix = f"[{suffix};{inserted_unit[:inserted_remainder]}]"

            if trimmed.start == trimmed.end:
                return var_dict("insertion", trimmed.start, trimmed.end, suffix)
                # return f"{to_hgvs_position(trimmed.start, trimmed.end)}ins{suffix}"
            return var_dict("deletion_insertion", trimmed.start, trimmed.end, suffix)
            # return f"{to_hgvs_position(trimmed.start, trimmed.end)}delins{suffix}"

        # All other variants
        return other(trimmed)
        # return trimmed.to_hgvs(ref_seq)

    if not variants:
        return []

    if len(variants) == 1:
        return [hgvs(variants[0])]

    return [hgvs(variant) for variant in variants]


def _add_minimal(graph, reference, output, prefix=""):
    minimal_descriptions = []
    minimal_length = 100

    for variants in itertools.islice(graph.paths(), minimal_length):
        reference_variants = []
        for variant in variants:
            reference_variants.append(
                Variant(variant.start, variant.end, variant.sequence)
            )
        minimal_descriptions.append(
            f"{prefix}{to_hgvs_experimental(reference_variants, reference)}"
        )
    output["minimal_descriptions"] = minimal_descriptions
    if len(minimal_descriptions) == minimal_length:
        output["first_minimal"] = minimal_length


def _no_protein_support():
    return {
        "errors": [
            {
                "code": "ENOPROTEINSUPPORT",
                "details": "Protein descriptions not supported in this experimental service.",
            }
        ]
    }


def _algebra_variants(variants_delins, sequences):
    variants_algebra = []
    for variant in variants_delins:
        variants_algebra.append(
            Variant(get_start(variant), get_end(variant), get_inserted_sequence(variant, sequences))
        )
    return variants_algebra


def _only_variants(d, algebra_hgvs, supremal, local_supremals, ref_seq, graph):
    d.normalized_description = algebra_hgvs
    d.de_hgvs_model = {"variants": to_model(algebra_hgvs, "variants")}
    output = d.output()
    output["supremal"] = {
        "hgvs": supremal.to_hgvs(),
        "spdi": supremal.to_spdi(ref_seq),
    }
    output["view_corrected"] = view_delins(
        d.delins_model["variants"], d.corrected_model["variants"], d.get_sequences()
    )
    d_n = Description(
        description=d.normalized_description,
        only_variants=True,
        sequence=ref_seq,
    )
    d_n.to_delins()
    output["view_normalized"] = view_delins(
        d_n.delins_model["variants"], d.de_hgvs_model["variants"], d.get_sequences()
    )
    output["influence"] = {"min_pos": supremal.start, "max_pos": supremal.end}
    output["dot"] = "\n".join(to_dot(ref_seq, graph))
    _add_minimal(graph, ref_seq, output)
    return output


def _descriptions(d, algebra_hgvs, supremal, graph):
    ref_seq = d.references["reference"]["sequence"]["seq"]

    algebra_model = {
        "type": d.corrected_model["type"],
        "reference": {"id": d.corrected_model["reference"]["id"]},
        "coordinate_system": "g",
        "variants": to_model(algebra_hgvs, "variants"),
    }
    internal = to_internal_indexing(to_internal_coordinates(algebra_model, d.get_sequences()))
    if d.corrected_model.get("predicted"):
        internal["predicted"] = True
    delins = to_delins(internal)

    d.de_hgvs_internal_indexing_model = internal
    d.construct_de_hgvs_internal_indexing_model()
    d.construct_de_hgvs_coordinates_model()
    d.construct_normalized_description()
    d.construct_genomic_equivalent()
    # d.construct_equivalent()

    output = d.output()

    output["algebra"] = algebra_hgvs
    output["supremal"] = {
        "hgvs": f"{d.corrected_model['reference']['id']}:g.{supremal.to_hgvs()}",
        "spdi": supremal.to_spdi(d.corrected_model["reference"]["id"]),
    }

    output["view_corrected"] = view_delins(
        d.delins_model["variants"],
        d.corrected_model["variants"] if not d.is_inverted() else list(reversed(d.corrected_model["variants"])),
        d.get_sequences(),
        invert=d.is_inverted(),
    )
    output["view_normalized"] = view_delins(
        delins["variants"],
        d.de_hgvs_model["variants"] if not d.is_inverted() else list(reversed(d.de_hgvs_model["variants"])),
        d.get_sequences(),
        invert=d.is_inverted(),
    )
    output["dot"] = "\n".join(to_dot(ref_seq, graph))
    _add_minimal(graph, ref_seq, output, f"{d.corrected_model['reference']['id']}:g.")

    return output


def view_algebra_variants(variants, ref_seq, names=None):
    if names is None:
        names = [algebra_variant_to_name_model(v) for v in variants]
    return view_delins(
        [algebra_variant_to_delins(v) for v in variants],
        names,
        {"reference": ref_seq},
    )


def _normalize_alt(d_m):
    d = Description(description_model=d_m)
    d.to_delins()

    if d.corrected_model.get("type") == "description_protein":
        return _no_protein_support()

    if d.errors:
        return d.output()

    if d.only_equals() or d.no_operation():
        d.normalize_only_equals_or_no_operation()
        d.remove_superfluous_selector()
        return d.output()

    algebra_variants = _algebra_variants(d.delins_model["variants"], d.get_sequences())
    ref_seq = d.references["reference"]["sequence"]["seq"]

    algebra_extracted_variants, graph = extract_variants(ref_seq, algebra_variants)
    supremal = graph.supremal

    algebra_hgvs = to_hgvs_experimental(algebra_extracted_variants, ref_seq)
    # assert variants_to_description(to_hgvs_dict(algebra_extracted_variants, ref_seq)) != algebra_hgvs

    return _descriptions(d, algebra_hgvs, supremal, graph)


def normalize_alt(description, only_variants=False, sequence=None):
    def _protein():
        protein_selector_model = get_protein_selector_model(
            d.references["reference"]["annotations"],
            get_selector_id(d.de_hgvs_model),
        )
        if protein_selector_model:
            if has_recoding_translation_exception(protein_selector_model):
                output.setdefault("infos", []).append(
                    infos.in_frame_stop_codon(get_selector_id(d.de_hgvs_model))
                )
            p_d = get_protein_description(
                variants_to_delins(d.de_hgvs_internal_indexing_model["variants"]),
                d.references,
                protein_selector_model
            )
            protein = dict(
                zip(["description", "reference", "predicted"], p_d[:3]))
            if len(p_d) == 6:
                protein["position_first"] = p_d[3]
                protein["position_last_original"] = p_d[4]
                protein["position_last_predicted"] = p_d[5]
            output["protein"] = protein

    d = Description(description=description, only_variants=only_variants, sequence=sequence)
    d.to_delins()

    if d.corrected_model.get("type") == "description_protein":
        return _no_protein_support()

    if d.errors:
        return d.output()

    if d.only_equals() or d.no_operation():
        d.normalize_only_equals_or_no_operation()
        d.remove_superfluous_selector()
        return d.output()

    algebra_variants = _algebra_variants(d.delins_model["variants"], d.get_sequences())
    ref_seq = d.references["reference"]["sequence"]["seq"]

    algebra_extracted_variants, graph = extract_variants(ref_seq, algebra_variants)
    supremal = graph.supremal

    algebra_hgvs = to_hgvs_experimental(algebra_extracted_variants, ref_seq)
    # assert variants_to_description(to_hgvs_dict(algebra_extracted_variants, ref_seq)) == algebra_hgvs
    # assert to_model(algebra_hgvs, "variants") == to_hgvs_dict(algebra_extracted_variants, ref_seq)

    local_supremals = local_supremal(ref_seq, graph)
    if only_variants:
        output = _only_variants(d, algebra_hgvs, supremal, local_supremals, ref_seq, graph)
    else:
        output = _descriptions(d, algebra_hgvs, supremal, graph)

    if d.is_inverted():
        ref_seq_complement = reverse_complement(ref_seq)
        algebra_extracted_variants_complement = []
        len_seq = len(ref_seq)
        supremal_complement = Variant(len_seq - supremal.end, len_seq - supremal.start, reverse_complement(supremal.sequence))
        algebra_extracted_variants_complement.append(supremal_complement)
        algebra_extracted_variants_complement, graph_complement = extract_variants(ref_seq_complement, algebra_extracted_variants_complement)
        output["dot_complement"] = "\n".join(to_dot(ref_seq_complement, graph_complement))
        algebra_hgvs_complement = to_hgvs_experimental(algebra_extracted_variants_complement, ref_seq_complement)
        # assert to_model(algebra_hgvs_complement, "variants") == to_hgvs_dict(algebra_extracted_variants_complement, ref_seq_complement)
        algebra_model = {
            "type": d.corrected_model["type"],
            "reference": {"id": d.corrected_model["reference"]["id"]},
            "coordinate_system": "g",
            "variants": to_model(algebra_hgvs_complement, "variants"),
        }
        internal = to_internal_coordinates(algebra_model, {"reference": ref_seq_complement})
        selector_model = d.get_selector_model()
        coordinate_system = d.corrected_model.get("coordinate_system")
        new_exon = []
        new_cds = []
        for exon in reversed(selector_model.get("exon", [])):
            new_exon.append((len_seq - exon[1], len_seq - exon[0]))
        for cds in reversed(selector_model.get("cds", [])):
            new_cds.append((len_seq - cds[1], len_seq - cds[0]))
        selector_model["exon"] = new_exon
        selector_model["cds"] = new_cds
        selector_model["inverted"] = False
        crossmap = crossmap_to_hgvs_setup(coordinate_system, selector_model, True)
        hgvs_model = initialize_hgvs_model(internal, coordinate_system, d.get_selector_id())
        for point, path in yield_point_locations_for_main_reference(internal):
            set_by_path(hgvs_model, path, point_to_hgvs(point, **crossmap))
        if d.corrected_model.get("predicted"):
            hgvs_model["predicted"] = True
        if output.get("normalized_description"):
            output["normalized_description"] = model_to_string(hgvs_model)

    output["view_local_supremal"] = view_algebra_variants(local_supremals, ref_seq)

    output["influence"] = [(v.start, v.end) for v in local_supremals]

    if d.de_hgvs_model.get("coordinate_system") in ["c", "n"]:
        output["rna"] = dna_to_rna(description)
        if output["rna"].get("errors") is None:
            _protein()
    elif d.de_hgvs_model.get("coordinate_system") in ["r"]:
        m = rna_to_dna(description)
        if not m.get("errors"):
            _predicted_dna = _normalize_alt(m)
            if _predicted_dna.get("normalized_model"):
                _predicted_dna["normalized_model"]["predicted"] = True
                output["dna"] = {"description": model_to_string(_predicted_dna["normalized_model"])}
        else:
            output["dna"] = m
        _protein()
    return output


def normalize(description, only_variants=False, sequence=None):
    """
    Obtain the normalized (canonical) variant representation.
    """
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )
    d.normalize()
    output = d.output()
    return output


def delins_model(description, only_variants=False, sequence=None):
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )
    d.to_delins()
    output = d.output()
    if d.delins_model:
        for variant in d.delins_model["variants"]:
            if variant.get("inserted"):
                for inserted in variant.get("inserted"):
                    if not inserted.get("sequence"):
                        inserted["sequence"] = construct_sequence(
                            [inserted], d.get_sequences()
                        )
        output["delins_model"] = d.delins_model
    return output
