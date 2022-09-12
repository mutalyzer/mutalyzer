from algebra import Variant
from algebra.extractor import extract_supremal, to_hgvs
from algebra.lcs.all_lcs import edit, lcs_graph, traversal
from algebra.relations.supremal_based import find_supremal, spanning_variant
from algebra.utils import to_dot
from algebra.variants import patch
from mutalyzer_hgvs_parser import to_model

from mutalyzer.util import get_end, get_inserted_sequence, get_start

from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .description import Description
from .viewer import view_delins


def _add_dot(supremal, reference, output):
    ref_seq = reference[supremal.start : supremal.end]
    obs_seq = supremal.sequence
    _, lcs_nodes = edit(ref_seq, obs_seq)
    if len(lcs_nodes) < 50:
        root, _ = lcs_graph(ref_seq, obs_seq, lcs_nodes)
        output["dot"] = to_dot(reference, root)
        minimal_descriptions = []
        for variants in traversal(root):
            reference_variants = []
            for variant in variants:
                reference_variants.append(
                    Variant(
                        supremal.start + variant.start,
                        supremal.start + variant.end,
                        variant.sequence,
                    )
                )
            minimal_descriptions.append(to_hgvs(reference_variants, reference))
        if minimal_descriptions and len(minimal_descriptions) <= 100:
            output["minimal_descriptions"] = minimal_descriptions


def _no_protein_support():
    return {
        "errors": [
            {
                "code": "ENOPROTEINSUPPORT",
                "details": f"Protein descriptions not supported in this experimental service.",
            }
        ]
    }


def _algebra_variants(variants_delins, sequences):
    variants_algebra = []
    for variant in variants_delins:
        variants_algebra.append(
            Variant(
                get_start(variant),
                get_end(variant),
                get_inserted_sequence(variant, sequences),
            )
        )
    return variants_algebra


def _only_variants(d, algebra_hgvs, supremal, reference):
    d.normalized_description = algebra_hgvs
    d.de_hgvs_model = {"variants": to_model(algebra_hgvs, "variants")}
    output = d.output()
    output["supremal"] = {
        "hgvs": f"{reference}:g.{supremal.to_hgvs()}",
        "spdi": supremal.to_spdi(reference),
    }
    output["view_corrected"] = {
        "views": view_delins(
            d.delins_model["variants"], d.corrected_model["variants"], reference
        ),
        "seq_length": len(reference),
    }
    d_n = Description(
        description=d.normalized_description,
        only_variants=True,
        sequence=reference,
    )
    d_n.to_delins()
    output["view_normalized"] = {
        "views": view_delins(
            d_n.delins_model["variants"], d.de_hgvs_model["variants"], reference
        ),
        "seq_length": len(reference),
    }
    output["influence"] = {"min_pos": supremal.start, "max_pos": supremal.end}
    return output


def _descriptions(d, algebra_hgvs, supremal, reference):
    algebra_model = {
        "type": d.corrected_model["type"],
        "reference": {"id": d.corrected_model["reference"]["id"]},
        "coordinate_system": "g",
        "variants": to_model(algebra_hgvs, "variants"),
    }

    d.de_hgvs_internal_indexing_model = to_internal_indexing(
        to_internal_coordinates(algebra_model, d.get_sequences())
    )
    d.construct_de_hgvs_internal_indexing_model()
    d.construct_de_hgvs_coordinates_model()
    d.construct_normalized_description()

    output = d.output()
    output["algebra"] = algebra_hgvs
    output["supremal"] = {
        "hgvs": f"{d.corrected_model['reference']['id']}:g.{supremal.to_hgvs()}",
        "spdi": supremal.to_spdi(d.corrected_model["reference"]["id"]),
    }

    output["view_corrected"] = {
        "views": view_delins(
            d.delins_model["variants"], d.corrected_model["variants"], reference
        ),
        "seq_length": len(reference),
    }
    d_n = Description(d.normalized_description)
    d_n.to_delins()
    output["view_normalized"] = {
        "views": view_delins(
            d_n.delins_model["variants"], d.de_hgvs_model["variants"], reference
        ),
        "seq_length": len(reference),
    }
    return output


def name_check_alt(description, only_variants=False, sequence=None):

    # TODO: reverse strand shift?

    d = Description(
        description=description, only_variants=only_variants, sequence=sequence
    )

    d.to_delins()

    if d.errors:
        return d.output()

    if not only_variants and d.corrected_model["type"] == "description_protein":
        _no_protein_support()

    algebra_variants = _algebra_variants(d.delins_model["variants"], d.get_sequences())

    ref_seq = d.references["reference"]["sequence"]["seq"]
    obs_seq = patch(ref_seq, algebra_variants)
    supremal = find_supremal(
        ref_seq, spanning_variant(ref_seq, obs_seq, algebra_variants)
    )
    canonical = list(extract_supremal(ref_seq, supremal))

    algebra_hgvs = to_hgvs(canonical, ref_seq)

    if only_variants:
        output = _only_variants(d, algebra_hgvs, supremal, ref_seq)
    else:
        output = _descriptions(d, algebra_hgvs, supremal, ref_seq)

    output["influence"] = {"min_pos": supremal.start, "max_pos": supremal.end}

    _add_dot(supremal, ref_seq, output)

    return output


def name_check(description, only_variants=False, sequence=None):
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )

    d.normalize()

    return d.output()
