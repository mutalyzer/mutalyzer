"""Interfaces to obtain the normalized (canonical) variant representations or
the delins model of an input description."""

import itertools

from algebra import Variant
from algebra.extractor import extract as extract_variants
from algebra.extractor import local_supremal
from algebra.extractor import to_hgvs as to_hgvs_experimental
from algebra.utils import to_dot
from mutalyzer_hgvs_parser import to_model

from mutalyzer.util import get_inserted_sequence, get_location_length

from .algebra import algebra_variant_to_delins, algebra_variant_to_name_model
from .converter.to_delins import to_delins
from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .description import Description
from .util import construct_sequence, get_end, get_start, roll
from .viewer import view_delins
from .rna import dna_to_rna, rna_to_dna
from .protein import get_protein_description
from .reference import get_protein_selector_model
from .description_model import get_selector_id, model_to_string
from .converter.to_delins import to_delins, variants_to_delins


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


def _add_shift(internal, delins, reference):
    for i, v in enumerate(delins["variants"]):
        inserted_sequence = get_inserted_sequence(v, {"reference": reference})
        shift5 = 0
        if get_location_length(v) and not inserted_sequence:
            shift5, _ = roll(
                reference,
                get_start(v) + 1,
                get_end(v),
            )
        elif not get_location_length(v) and inserted_sequence:
            rolled_sequence = (
                reference[: get_start(v)] + inserted_sequence + reference[get_end(v) :]
            )
            shift5, _ = roll(
                rolled_sequence,
                get_start(v) + 1,
                get_end(v) + len(inserted_sequence),
            )
        internal["variants"][i]["location"]["start"]["shift"] = shift5
        internal["variants"][i]["location"]["end"]["shift"] = shift5


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
    _add_shift(internal, delins, ref_seq)

    d.de_hgvs_internal_indexing_model = internal
    d.construct_de_hgvs_internal_indexing_model()
    d.construct_de_hgvs_coordinates_model()
    d.construct_normalized_description()
    d.construct_equivalent()

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

    return _descriptions(d, algebra_hgvs, supremal, graph)


def normalize_alt(description, only_variants=False, sequence=None):
    def _protein():
        protein_selector_model = get_protein_selector_model(
            d.references["reference"]["annotations"],
            get_selector_id(d.de_hgvs_model),
        )
        if protein_selector_model:
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

    local_supremals = local_supremal(ref_seq, graph)
    if only_variants:
        output = _only_variants(d, algebra_hgvs, supremal, local_supremals, ref_seq, graph)
    else:
        output = _descriptions(d, algebra_hgvs, supremal, graph)

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
    d.get_chromosomal_descriptions()
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
