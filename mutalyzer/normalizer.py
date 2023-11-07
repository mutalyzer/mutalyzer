"""Interfaces to obtain the normalized (canonical) variant representations or
the delins model of an input description."""

import bisect
import itertools
from collections import deque

from algebra import Variant
from algebra.extractor import extract as extract_variants
from algebra.extractor import local_supremal
from algebra.extractor import to_hgvs as to_hgvs_experimental
from algebra.lcs.all_lcs import dfs_traversal
from algebra.utils import to_dot
from algebra.variants import patch
from mutalyzer_crossmapper import NonCoding
from mutalyzer_hgvs_parser import to_model

from mutalyzer.util import get_inserted_sequence, get_location_length

from .algebra import algebra_variant_to_delins, algebra_variant_to_name_model
from .converter.to_delins import to_delins
from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
from .description import Description
from .util import construct_sequence, get_end, get_start, roll
from .viewer import view_delins


def add_pre_edges(root):
    root.pre_edges = []
    visited = {root}
    queue = deque([root])
    while queue:
        node = queue.popleft()
        if not node.edges:
            sink = node
        for succ, variant in node.edges:
            if succ not in visited:
                succ.pre_edges = []
                visited.add(succ)
                queue.append(succ)
            succ.pre_edges.append((node, variant))
    return sink


def get_sides_limits(root, paths=True):
    sink = add_pre_edges(root)

    right = max(root.edges, key=lambda x: x[1].start)
    right_pos = right[1].start

    left = min(sink.pre_edges, key=lambda x: x[1].end)
    left_pos = left[1].end

    right_path = []
    left_path = []

    if paths:
        right_path.append(right[1])
        node = right[0]
        while node.edges:
            right_path.append(node.edges[0][1])
            node = node.edges[0][0]

        left_path.append(left[1])
        node = left[0]
        while node.pre_edges:
            left_path.append(node.pre_edges[0][1])
            node = node.pre_edges[0][0]

    return (right_pos, right_path), (left_pos, left_path)


def _to_dot(reference, root, offset=0):
    """The LCS graph in Graphviz DOT format."""

    def traverse():
        # breadth-first traversal
        node_count = 0
        visited = {root: node_count}
        queue = deque([root])
        while queue:
            node = queue.popleft()
            if not node.edges:
                yield f'"s{visited[node]}"' + "[shape=doublecircle]"
            for succ, variant in node.edges:
                if succ not in visited:
                    node_count += 1
                    visited[succ] = node_count
                    queue.append(succ)
                hgvs_variant = Variant(
                    variant.start + offset,
                    variant.end + offset,
                    variant.sequence,
                ).to_hgvs(reference)
                yield (
                    f'"s{visited[node]}" -> "s{visited[succ]}"'
                    f' [label="{hgvs_variant}"];'
                )

    return (
            "digraph {\n"
            "    rankdir=LR\n"
            "    edge[fontname=monospace]\n"
            "    node[shape=circle]\n"
            "    si[shape=point]\n"
            "    si->" + f'"s{0}"' + "\n"
                                     "    " + "\n    ".join(traverse()) + "\n}"
    )


def get_rna_limits(root, ref_seq, offset=0):
    left, right = get_sides_limits(root)
    print(left, right)

    left_variants = [Variant(v.start + offset, v.end + offset, v.sequence) for v in left[1]]
    right_variants = [Variant(v.start + offset, v.end + offset, v.sequence) for v in right[1]]

    left_model = to_model("[" + ";".join([v.to_hgvs(ref_seq) for v in left_variants]) + "]", "variants")
    right_model = to_model("[" + ";".join([v.to_hgvs(ref_seq) for v in right_variants]) + "]", "variants")
    # print(str(left[0] + offset) + ": [" + ";".join([v.to_hgvs(ref_seq) for v in left_variants]) + "]")
    # print(str(right[0] + offset) + ": [" + ";".join([v.to_hgvs(ref_seq) for v in right_variants]) + "]")
    return (left[0] + offset, left_model), (right[0] + offset, right_model)


def post_dominators(node, start, visited, reference):
    if node in visited:
        visited[node]["start"] = max(visited[node]["start"], start)
        return visited

    visited[node] = {
        "post": {node},
        "start": start,
        "end": len(reference),
    }

    post = set()
    for child, variant in node.edges:
        post_dominators(child, variant.end, visited)
        if not post:
            post = visited[child]["post"]
        post = post.intersection(visited[child]["post"])

        visited[node]["end"] = min(visited[node]["end"], variant.start)

    visited[node]["post"] = post.union(visited[node]["post"])
    return visited


def _add_minimal(root, reference, output, prefix=""):
    minimal_descriptions = []
    minimal_length = 100
    for variants in itertools.islice(dfs_traversal(root), minimal_length):
        reference_variants = []
        for variant in variants:
            reference_variants.append(Variant(variant.start, variant.end, variant.sequence))
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
            Variant(
                get_start(variant),
                get_end(variant),
                get_inserted_sequence(variant, sequences),
            )
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


def _only_variants(d, algebra_hgvs, supremal, ref_seq, root):
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
    output["dot"] = "\n".join(to_dot(ref_seq, root))
    _add_minimal(root, ref_seq, output)
    return output


def get_position_type(position, exons, len_ss=2):
    def _output_exon(index):
        exon_index = index // 2
        left_shift = position - exons[exon_index][0] + 1
        right_shift = position - exons[exon_index][1]
        if abs(left_shift) <= len_ss:
            return index, left_shift
        if abs(right_shift) <= len_ss:
            return index, right_shift
        return index, 0

    def _output_intron(index, offset):
        if abs(offset) <= len_ss:
            return index, offset
        else:
            return index, 0

    x = NonCoding(exons).coordinate_to_noncoding
    flattened_exons = [e for exon in exons for e in exon]
    position_x = x(position)

    if position_x[1] == 0:
        return _output_exon(bisect.bisect_right(flattened_exons, position))
    elif 0 < abs(position_x[1]) <= len_ss and position_x[1] > 0:
        return _output_intron(bisect.bisect_right(flattened_exons, position), position_x[1])
    else:
        return _output_intron(bisect.bisect_left(flattened_exons, position), position_x[1])


def construct_rna_description(d, ref_seq, root, algebra_variants):
    # NG_012337.3(NM_003002.4):c.172_175dup
    # NG_012337.3(NM_003002.4):c.[310_314+4dup]
    # NG_012337.3(NM_00-3002.4):c.[300del;310_314+4dup]
    exons = d.get_selector_model()["exon"]
    print("=====")
    print("exons:")
    print(exons)

    exon_margin = 2
    intron_margin = 4

    local_sup = local_supremal(ref_seq, patch(ref_seq, algebra_variants), root)
    print("local supremals:")
    for sup in local_sup:
        print("----")
        print(ref_seq[sup.start:sup.end])
        print(sup)
        sup_start_index, sup_start_offset = get_position_type(sup.start, exons, exon_margin)
        sup_end_index, sup_end_offset = get_position_type(sup.end, exons, intron_margin)
        print(sup_start_index, sup_end_index)
        if sup_end_index - sup_start_index == 1:
            print("we have to check this")
            print("rna limits")
            right_push, left_push = get_sides_limits(root)
            print(right_push)
            print(left_push)
            if sup_start_index % 2 == 0:
                print("intron - exon")
                print(exons[sup_start_index//2][1])
                if left_push[0] > exons[sup_start_index//2][0] - intron_margin:
                    print("it can be pushed into the intron")
                if right_push[0] > exons[sup_start_index//2][0] + exon_margin:
                    print("it can be pushed into the exon")
            else:
                print("exon - intron")
                print(exons[sup_start_index//2][1])
                if right_push[0] > exons[sup_start_index//2][1] + intron_margin - 1:
                    print("it can be pushed into the intron")
                if left_push[0] > exons[sup_start_index//2][1] - exon_margin:
                    print("it can be pushed into the exon")


def _descriptions(d, algebra_variants, algebra_hgvs, supremal, ref_seq, root):
    algebra_model = {
        "type": d.corrected_model["type"],
        "reference": {"id": d.corrected_model["reference"]["id"]},
        "coordinate_system": "g",
        "variants": to_model(algebra_hgvs, "variants"),
    }
    internal = to_internal_indexing(to_internal_coordinates(algebra_model, d.get_sequences()))
    delins = to_delins(internal)
    _add_shift(internal, delins, ref_seq)

    d.de_hgvs_internal_indexing_model = internal
    d.construct_de_hgvs_internal_indexing_model()
    d.construct_de_hgvs_coordinates_model()
    d.construct_normalized_description()
    d.construct_equivalent()

    # if d.de_hgvs_model.get("coordinate_system") in ["c", "n"]:
    #     construct_rna_description(d, ref_seq, root, algebra_variants)

    d.construct_protein_description()

    output = d.output()
    output["algebra"] = algebra_hgvs
    output["supremal"] = {
        "hgvs": f"{d.corrected_model['reference']['id']}:g.{supremal.to_hgvs()}",
        "spdi": supremal.to_spdi(d.corrected_model["reference"]["id"]),
    }

    output["view_corrected"] = view_delins(
        d.delins_model["variants"],
        d.corrected_model["variants"],
        d.get_sequences(),
        invert=d.is_inverted(),
    )
    output["view_normalized"] = view_delins(
        delins["variants"],
        d.de_hgvs_model["variants"],
        d.get_sequences(),
        invert=d.is_inverted(),
    )
    output["dot"] = "\n".join(to_dot(ref_seq, root))
    _add_minimal(root, ref_seq, output, f"{d.corrected_model['reference']['id']}:g.")
    return output


def normalize_alt(description, only_variants=False, sequence=None):
    d = Description(
        description=description, only_variants=only_variants, sequence=sequence
    )
    d.to_delins()
    if d.corrected_model.get("type") == "description_protein":
        p_d = Description(
            description=description, only_variants=only_variants, sequence=sequence
        )
        p_d.normalize()
        return p_d.output()

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

    algebra_extracted_variants, supremal, root = extract_variants(ref_seq, algebra_variants)

    algebra_hgvs = to_hgvs_experimental(algebra_extracted_variants, ref_seq)

    if only_variants:
        output = _only_variants(d, algebra_hgvs, supremal, ref_seq, root)
    else:
        output = _descriptions(d, algebra_extracted_variants, algebra_hgvs, supremal, ref_seq, root)

    local_supremals = local_supremal(ref_seq, patch(ref_seq, algebra_extracted_variants), root)

    output["view_local_supremal"] = view_delins(
        [algebra_variant_to_delins(v) for v in local_supremals],
        [algebra_variant_to_name_model(v) for v in local_supremals],
        {"reference": ref_seq},
    )

    output["influence"] = {"min_pos": supremal.start, "max_pos": supremal.end}

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
