import bisect
import json
from collections import deque

from algebra import LCSgraph, Variant
from algebra.extractor import extract as extract_variants
from algebra.extractor import local_supremal as get_local_supremal
from algebra.variants import patch, to_hgvs
from Bio.Seq import Seq
from mutalyzer_crossmapper import NonCoding
from mutalyzer_hgvs_parser import to_model

from .algebra import (
    algebra_variant_to_delins,
    delins_to_algebra_variant,
    delins_to_algebra_variants
)
from .converter.to_hgvs_coordinates import to_hgvs_locations
from .converter.to_rna import to_rna_reference_model, to_rna_variants
from .converter.variants_de_to_hgvs import de_to_hgvs
from .description import Description
from .description_model import (
    get_reference_id,
    variant_to_description,
    variants_to_description,
    model_to_string
)


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
        return _output_intron(
            bisect.bisect_right(flattened_exons, position), position_x[1]
        )
    else:
        return _output_intron(
            bisect.bisect_left(flattened_exons, position), position_x[1]
        )


def algebra_variants_to_hgvs(algebra_variants):
    variants = []
    for v in sorted(algebra_variants):
        variant = {
            "location": {
                "type": "range",
                "start": {"type": "point", "position": v.start},
                "end": {"type": "point", "position": v.end},
            }
        }
        if v.start == v.end:
            variant["type"] = "insertion"
            variant["inserted"] = [{"sequence": v.sequence, "source": "description"}]
        elif not v.sequence:
            variant["type"] = "deletion"
            variant["inserted"] = []
        else:
            variant["type"] = "deletion_insertion"
            variant["inserted"] = [{"sequence": v.sequence, "source": "description"}]
        variants.append(variant)
    return variants


def extracted_to_hgvs_selector(variants, d, to_selector_id):
    extracted_model = {
        "reference": {"id": d.corrected_model["reference"]["id"]},
        "variants": algebra_variants_to_hgvs(variants),
    }
    de_hgvs_model = to_hgvs_locations(
        model=extracted_model,
        references=d.references,
        to_selector_id=to_selector_id,
        degenerate=True,
    )
    return de_hgvs_model


def get_rna_sequences(reference_model, selector_id):

    rna_reference_model = to_rna_reference_model(reference_model, selector_id)
    return {
        reference_model["id"]: rna_reference_model,
        "reference": rna_reference_model,
    }


def get_rna_sequences(d):
    rna_references = get_rna_reference_models(d)
    return {
        k: str(Seq(rna_references[k]["sequence"]["seq"]).transcribe().lower())
        for k in rna_references
    }


def get_rna_reference_models(d):
    reference_model = d.references["reference"]
    reference_id = get_reference_id(d.corrected_model)
    selector_id = d.get_selector_id()

    rna_reference_model = to_rna_reference_model(reference_model, selector_id)
    return {reference_id: rna_reference_model, "reference": rna_reference_model}


def get_rna_variants(d, variants):
    delins_variants = [algebra_variant_to_delins(v) for v in variants]
    rna_delins = to_rna_variants(
        delins_variants,
        d.get_sequences(),
        d.get_selector_model(),
    )
    rna_reference_models = get_rna_reference_models(d)
    rna_ref_seq = rna_reference_models["reference"]["sequence"]["seq"]
    rna_algebra_variants = [
        delins_to_algebra_variant(v, get_rna_sequences(d)) for v in rna_delins
    ]
    rna_algebra_extracted_variants, *_ = extract_variants(
        rna_ref_seq, rna_algebra_variants
    )
    rna_variants_coordinate = de_to_hgvs(
        [algebra_variant_to_delins(v) for v in rna_algebra_extracted_variants],
        {k: rna_reference_models[k]["sequence"]["seq"] for k in rna_reference_models},
    )
    rna_model = to_hgvs_locations(
        {
            "reference": d.corrected_model["reference"],
            "coordinate_system": "i",
            "variants": rna_variants_coordinate,
        },
        rna_reference_models,
        d.de_hgvs_model.get("coordinate_system"),
        d.get_selector_id(),
        True,
    )
    rna_model["coordinate_system"] = "r"
    rna_model["predicted"] = True
    rna_algebra_hgvs = [variant_to_description(v) for v in rna_model["variants"]]
    return rna_algebra_hgvs


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


def get_rna_limits(graph, ref_seq, offset=0):
    left, right = get_sides_limits(graph)
    # print(left, right)

    left_variants = [
        Variant(v.start + offset, v.end + offset, v.sequence) for v in left[1]
    ]
    right_variants = [
        Variant(v.start + offset, v.end + offset, v.sequence) for v in right[1]
    ]

    left_model = to_model(
        "[" + ";".join([v.to_hgvs(ref_seq) for v in left_variants]) + "]", "variants"
    )
    right_model = to_model(
        "[" + ";".join([v.to_hgvs(ref_seq) for v in right_variants]) + "]", "variants"
    )
    return (left[0] + offset, left_model), (right[0] + offset, right_model)


def get_sides_limits(graph, paths=True):
    root = graph._source
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


def _genomic_and_coding(algebra_variants, d, selector_id):
    ref_seq = d.references["reference"]["sequence"]["seq"]
    return {
        "genomic": to_hgvs(algebra_variants, ref_seq),
        "coding": [
            variant_to_description(v)
            for v in extracted_to_hgvs_selector(algebra_variants, d, selector_id)[
                "variants"
            ]
        ],
    }


def predict_rna(d, local_supremals):
    # NG_012337.3(NM_003002.4):c.172_175dup
    # NG_012337.3(NM_003002.4):c.[310_314+4dup]
    # NG_012337.3(NM_003002.4):c.[300del;310_314+4dup]

    selector_id = d.get_selector_id()
    ref_seq = d.references["reference"]["sequence"]["seq"]
    exons = d.get_selector_model()["exon"]

    exon_margin = 2
    intron_margin = 4

    status = {
        "exon_margin": exon_margin,
        "intron_margin": intron_margin,
        "local_supremals": {},
    }
    for i, sup in enumerate(local_supremals):
        _, local_root = extract_variants(ref_seq, [sup])
        sup_start_index, sup_start_offset = get_position_type(sup.start, exons, exon_margin)
        sup_end_index, sup_end_offset = get_position_type(sup.end, exons, intron_margin)
        print(sup_start_index, sup_start_offset)
        print(sup_end_index, sup_end_offset)
        splice_affected = False
        sup_status = {}
        if sup_end_index - sup_start_index == 1:
            splice_affected = True
            right_push, left_push = get_sides_limits(local_root)
            if sup_start_index % 2 == 0:
                # intron - exon
                sup_status["between"] = "intron - exon"
                if left_push[0] < exons[sup_start_index // 2][0] - intron_margin:
                    # it can be pushed into the intron
                    sup_status["push_intron"] = left_push[1]
                if right_push[0] > exons[sup_start_index // 2][0] + exon_margin:
                    # it can be pushed into the exon
                    # print(right_push)
                    sup_status["push_exon"] = right_push[1]
                    # print("\n\n\n\ssfd")
                    # print(sup_status["push_exon"])
            else:
                # exon - intron
                sup_status["between"] = "exon - intron"
                if right_push[0] > exons[sup_start_index // 2][1] + intron_margin - 1:
                    # it can be pushed into the intron
                    sup_status["push_intron"] = right_push[1]
                if left_push[0] < exons[sup_start_index // 2][1] - exon_margin:
                    # it can be pushed into the exon
                    sup_status["push_exon"] = left_push[1]
        # print(extract_variants(ref_seq, [sup]))
        # print(_genomic_and_coding(extract_variants(ref_seq, [sup])[0], d, selector_id))
        sup_status["hgvs"] = _genomic_and_coding(
            extract_variants(ref_seq, [sup])[0], d, selector_id
        )
        sup_status["splice_affected"] = splice_affected
        sup_status["supremal"] = _genomic_and_coding([sup], d, selector_id)
        status["local_supremals"][i] = sup_status

    print(json.dumps(status, indent=2))

    rna_description_possible = True
    for i, sup_status in status["local_supremals"].items():

        if sup_status.get("splice_affected"):
            if sup_status.get("push_intron") and sup_status.get("push_exon") is None:
                # Everything in intron, so no description at the RNA level.
                # print("Only an intronic variant.")
                sup_status["push_intron"] = _genomic_and_coding(
                    sup_status["push_intron"], d, selector_id
                )
            elif sup_status.get("push_intron") is None and sup_status.get("push_exon"):
                # print("see what can be done with the exon:", sup_status.get("push_exon"))
                sup_status["rna"] = get_rna_variants(d, sup_status.get("push_exon"))
                sup_status["push_exon"] = _genomic_and_coding(
                    sup_status["push_exon"], d, selector_id
                )
            elif sup_status.get("push_intron") and sup_status.get("push_exon"):
                sup_status["push_exon"] = _genomic_and_coding(
                    sup_status["push_exon"], d, selector_id
                )
                sup_status["push_intron"] = _genomic_and_coding(
                    sup_status["push_intron"], d, selector_id
                )
                rna_description_possible = False
            else:
                # The splice is always affected.
                # print("Nothing possible.")
                rna_description_possible = False
        else:
            sup_status["rna"] = get_rna_variants(d, [local_supremals[i]])

    if rna_description_possible:
        rna_description = []
        for i, sup_status in status["local_supremals"].items():
            if sup_status.get("rna"):
                rna_description.extend(sup_status["rna"])
        description = d.corrected_model["reference"]["id"]
        if d.get_selector_id():
            description += f"({d.get_selector_id()})"
        if len(rna_description) > 1:
            description += f":r.([{';'.join(rna_description)}])"
        else:
            description += f":r.({';'.join(rna_description)})"
        status["description"] = description

    # print(json.dumps(status, indent=2))
    return status


def _splice_sites_affected(exons, local_supremal, exon_margin=2, intron_margin=4):
    """
    Check if the local supremal variants touch the splice sites
    within the exon/intron margins.
    """
    for i, sup in enumerate(local_supremal):
        sup_start_index, sup_start_offset = get_position_type(sup.start, exons, exon_margin)
        sup_end_index, sup_end_offset = get_position_type(sup.end, exons, intron_margin)
        if sup_end_index - sup_start_index == 1:
            return True
    return False


def _to_rna_variants(variants, exons):
    x = NonCoding(exons).coordinate_to_noncoding
    rna_variants = []
    for variant in variants:
        start = x(variant.start)[0] - 1
        end = x(variant.end)[0] + x(variant.end)[1] - 1
        rna_variants.append(Variant(start, end, str(Seq(variant.sequence).transcribe().lower())))
    return rna_variants


def dna_to_rna(description):
    d = Description(description)
    d.to_delins()
    d.print_models_summary()
    if d.errors:
        return d.errors

    ref_seq = d.references["reference"]["sequence"]["seq"]
    algebra_variants = delins_to_algebra_variants(d.delins_model["variants"], d.get_sequences())

    algebra_extracted_variants, graph = extract_variants(ref_seq, algebra_variants)
    local_supremal = get_local_supremal(ref_seq, graph)

    if not _splice_sites_affected(d.get_selector_model()["exon"], local_supremal):
        print("we can do something")
        rna_sliced_variants = _to_rna_variants(algebra_extracted_variants, d.get_selector_model()["exon"])

        rna_reference_models = get_rna_reference_models(d)
        rna_ref_seq = rna_reference_models["reference"]["sequence"]["seq"]
        rna_algebra_extracted_variants, *_ = extract_variants(
            rna_ref_seq, rna_sliced_variants
        )
        rna_variants_coordinate = de_to_hgvs(
            [algebra_variant_to_delins(v) for v in
             rna_algebra_extracted_variants],
            {k: rna_reference_models[k]["sequence"]["seq"] for k in
             rna_reference_models},
        )

        rna_model = to_hgvs_locations(
            {
                "reference": d.corrected_model["reference"],
                "coordinate_system": "i",
                "variants": rna_variants_coordinate,
            },
            rna_reference_models,
            d.de_hgvs_model.get("coordinate_system"),
            d.get_selector_id(),
            True,
        )
        rna_model["coordinate_system"] = "r"
        rna_model["predicted"] = True
        print(rna_sliced_variants)
        print(model_to_string(rna_model))
        return model_to_string(rna_model)
    else:
        print("not really")
