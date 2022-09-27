import itertools
from collections import deque

from algebra import Variant
from algebra.extractor import to_hgvs
from algebra.extractor.extractor import canonical
from algebra.lcs.all_lcs import edit, lcs_graph, traversal
from algebra.relations.supremal_based import find_supremal, spanning_variant
from algebra.variants import patch
from mutalyzer_hgvs_parser import to_model

from mutalyzer.util import get_end, get_inserted_sequence, get_start, get_location_length

from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .converter.to_delins import to_delins
from .description import Description
from .viewer import view_delins
from .util import roll, get_start, get_end


def _to_dot(reference, root, offset):
    def traverse():
        visited = {root}
        queue = deque([root])
        while queue:
            node = queue.popleft()
            for succ, variant in node.edges:
                if variant:
                    hgvs_variant = [
                        Variant(
                            variant[0].start + offset,
                            variant[0].end + offset,
                            variant[0].sequence,
                        )
                    ]
                    style = ""
                else:
                    hgvs_variant = []
                    style = 'style = "dashed"'
                yield (
                    f'"{node.row}_{node.col}" [shape=circle, label=""];',
                    f'"{succ.row}_{succ.col}" [shape=circle, label=""];',
                    f' "{node.row}_{node.col}" -> "{succ.row}_{succ.col}"'
                    f' [label="{to_hgvs(hgvs_variant, reference)}" {style}];',
                )
                if succ not in visited:
                    visited.add(succ)
                    queue.append(succ)

    nodes = set()
    edges = []
    for node_start, node_end, edge in traverse():
        nodes.add(node_start)
        nodes.add(node_end)
        edges.append(edge)
        if len(nodes) > 100:
            return None
    elements = list(nodes) + edges
    return "digraph {\n    " + "\n    ".join(elements) + "\n}"


def _add_dot(supremal, root, reference, output, prefix=""):
    dot = _to_dot(reference, root, supremal.start)
    if dot:
        output["dot"] = _to_dot(reference, root, supremal.start)
    minimal_descriptions = []
    minimal_length = 100
    for variants in itertools.islice(traversal(root), minimal_length):
        reference_variants = []
        for variant in variants:
            reference_variants.append(
                Variant(
                    supremal.start + variant.start,
                    supremal.start + variant.end,
                    variant.sequence,
                )
            )
        minimal_descriptions.append(f"{prefix}{to_hgvs(reference_variants, reference)}")
    output["minimal_descriptions"] = minimal_descriptions
    if len(minimal_descriptions) == minimal_length:
        output["first_minimal"] = minimal_length


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


def _add_shift(internal, delins, reference):
    print("\n\n----")
    for i, v in enumerate(delins["variants"]):
        inserted_sequence = get_inserted_sequence(v, {"reference": reference})

        shift3 = 0
        shift5 = 0

        if get_location_length(v) and not inserted_sequence:
            shift5, shift3 = roll(
                reference,
                get_start(v) + 1,
                get_end(v),
            )
        elif not get_location_length(v) and inserted_sequence:
            rolled_sequence = (reference[: get_start(v)] + inserted_sequence + reference[get_end(v):])
            shift5, shift3 = roll(
                rolled_sequence,
                get_start(v) + 1,
                get_end(v) + len(inserted_sequence),
            )
        print(shift5, shift3)
        internal["variants"][i]["location"]["start"]["shift"] = shift5
        internal["variants"][i]["location"]["end"]["shift"] = shift5


def _only_variants(d, algebra_hgvs, supremal, ref_seq, root):
    d.normalized_description = algebra_hgvs
    d.de_hgvs_model = {"variants": to_model(algebra_hgvs, "variants")}
    output = d.output()
    output["supremal"] = {
        "hgvs": supremal.to_hgvs(),
        "spdi": supremal.to_spdi(ref_seq),
    }
    output["view_corrected"] = {
        "views": view_delins(
            d.delins_model["variants"], d.corrected_model["variants"],  d.get_sequences()
        ),
        "seq_length": len(ref_seq),
    }
    d_n = Description(
        description=d.normalized_description,
        only_variants=True,
        sequence=ref_seq,
    )
    d_n.to_delins()
    output["view_normalized"] = {
        "views": view_delins(
            d_n.delins_model["variants"], d.de_hgvs_model["variants"], d.get_sequences()
        ),
        "seq_length": len(ref_seq),
    }
    output["influence"] = {"min_pos": supremal.start, "max_pos": supremal.end}
    _add_dot(supremal, root, ref_seq, output)
    return output


def _descriptions(d, algebra_hgvs, supremal, ref_seq, root):
    algebra_model = {
        "type": d.corrected_model["type"],
        "reference": {"id": d.corrected_model["reference"]["id"]},
        "coordinate_system": "g",
        "variants": to_model(algebra_hgvs, "variants"),
    }
    print("\n--- descriptions")
    internal = to_internal_indexing(
        to_internal_coordinates(algebra_model, d.get_sequences())
    )
    delins = to_delins(internal)
    import json
    print(json.dumps(delins, indent=2))
    _add_shift(internal, delins, ref_seq)
    print(json.dumps(internal, indent=2))

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

    output["view_corrected"] = {
        "views": view_delins(
            d.delins_model["variants"], d.corrected_model["variants"],  d.get_sequences()
        ),
        "seq_length": len(ref_seq),
    }
    d_n = Description(d.normalized_description)
    d_n.to_delins()
    output["view_normalized"] = {
        "views": view_delins(
            d_n.delins_model["variants"], d.de_hgvs_model["variants"], d.get_sequences()
        ),
        "seq_length": len(ref_seq),
    }
    _add_dot(
        supremal, root, ref_seq, output, f"{d.corrected_model['reference']['id']}:g."
    )
    return output


def name_check_alt(description, only_variants=False, sequence=None):

    # TODO: reverse strand shift?

    d = Description(
        description=description, only_variants=only_variants, sequence=sequence
    )

    d.to_delins()

    if d.errors:
        return d.output()
    if d.only_equals() or d.no_operation():
        d.normalize_only_equals_or_no_operation()
        d.remove_superfluous_selector()
        return d.output()

    if not only_variants and d.corrected_model["type"] == "description_protein":
        _no_protein_support()

    algebra_variants = _algebra_variants(d.delins_model["variants"], d.get_sequences())

    ref_seq = d.references["reference"]["sequence"]["seq"]
    obs_seq = patch(ref_seq, algebra_variants)
    supremal = find_supremal(
        ref_seq, spanning_variant(ref_seq, obs_seq, algebra_variants)
    )

    supremal_ref_seq = ref_seq[supremal.start : supremal.end]
    supremal_obs_seq = supremal.sequence

    _, lcs_nodes = edit(supremal_ref_seq, supremal_obs_seq)
    root, _ = lcs_graph(supremal_ref_seq, supremal_obs_seq, lcs_nodes)

    algebra_hgvs = to_hgvs(
        list(
            [
                Variant(
                    supremal.start + variant.start,
                    supremal.start + variant.end,
                    variant.sequence,
                )
                for variant in canonical(supremal_obs_seq, root)
            ]
        ),
        ref_seq,
    )

    if only_variants:
        output = _only_variants(d, algebra_hgvs, supremal, ref_seq, root)

    else:
        output = _descriptions(d, algebra_hgvs, supremal, ref_seq, root)

    output["influence"] = {"min_pos": supremal.start, "max_pos": supremal.end}

    return output


def name_check(description, only_variants=False, sequence=None):
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )

    d.normalize()

    return d.output()
