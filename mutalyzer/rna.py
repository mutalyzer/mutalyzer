import bisect
import json
from collections import deque
from copy import deepcopy

from algebra import LCSgraph, Variant
from algebra.extractor import extract as extract_variants
from algebra.extractor import local_supremal as get_local_supremal
from algebra.extractor import to_hgvs
from Bio.Seq import Seq
from mutalyzer_crossmapper import Coding, Genomic, NonCoding
from mutalyzer_hgvs_parser import to_model
from mutalyzer_mutator.util import reverse_complement
from mutalyzer_retriever.reference import get_reference_mol_type
from mutalyzer_retriever.retriever import extract_feature_model

from .algebra import (
    algebra_variant_to_delins,
    delins_to_algebra,
    delins_to_algebra_variant,
)
from .converter.to_hgvs_coordinates import coding_to_point, to_hgvs_locations
from .converter.to_rna import to_rna_variants
from .converter.variants_de_to_hgvs import de_to_hgvs
from .description import Description
from .description_model import (
    get_reference_id,
    model_to_string,
    variant_to_description,
    yield_point_locations_for_main_reference,
    yield_ranges_main_reference,
)
from .errors import splice_site
from .reference import (
    get_coordinate_system_from_selector_id,
    get_internal_selector_model,
    slice_to_selector,
    yield_locations,
)
from .util import get_end, get_start, set_by_path, set_end, set_start


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

    return _output_intron(bisect.bisect_left(flattened_exons, position), position_x[1])


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

    rna_reference_model = to_rna_reference_model(reference_model, selector_id, d.is_inverted())
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
                    sup_status["push_exon"] = right_push[1]
            else:
                # exon - intron
                sup_status["between"] = "exon - intron"
                if right_push[0] > exons[sup_start_index // 2][1] + intron_margin - 1:
                    # it can be pushed into the intron
                    sup_status["push_intron"] = right_push[1]
                if left_push[0] < exons[sup_start_index // 2][1] - exon_margin:
                    # it can be pushed into the exon
                    sup_status["push_exon"] = left_push[1]
        sup_status["hgvs"] = _genomic_and_coding(
            extract_variants(ref_seq, [sup])[0], d, selector_id
        )
        sup_status["splice_affected"] = splice_affected
        sup_status["supremal"] = _genomic_and_coding([sup], d, selector_id)
        status["local_supremals"][i] = sup_status

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
        for _, sup_status in status["local_supremals"].items():
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

    return status


def _splice_sites_affected(exons, local_supremal, exon_margin=2, intron_margin=4):
    """
    Check if the local supremal variants touch the splice sites
    within the exon/intron margins.
    """
    for sup in local_supremal:
        sup_start_index, sup_start_offset = get_position_type(sup.start, exons, exon_margin)
        sup_end_index, sup_end_offset = get_position_type(sup.end, exons, intron_margin)
        if sup_end_index - sup_start_index == 1:
            return True
        elif (
            sup_end_index != sup_start_index
            and (sup_end_index - sup_start_index) % 2 == 0
            and sup.sequence
        ):
            return True
        for exon in exons:
            if exon[0] == sup.start == sup.end or exon[1] == sup.start == sup.end:
                return True
    return False


def _to_rna_variants(variants, exons):
    """
    Convert algebra variants locations to the rna (exons) slices.
    """
    x = NonCoding(exons).coordinate_to_noncoding
    rna_variants = []
    for variant in variants:

        start = variant.start
        end = variant.end

        intron_start_i = get_position_type(start, exons)[0]
        intron_end_i = get_position_type(end, exons)[0]

        if intron_start_i % 2 == 0:
            start = exons[intron_start_i // 2][0]
        if intron_end_i % 2 == 0:
            end = exons[intron_start_i // 2][1]

        start = x(start)[0] - 1
        end = x(end)[0] + x(end)[1] - 1

        rna_variants.append(Variant(start, end, str(Seq(variant.sequence).transcribe().lower())))

    return rna_variants


def to_rna_reference_model(reference_model, selector_id, inverted, transcribe=True):
    """
    Get the RNA reference model of the provided selector.

    1. Extract the tree corresponding to the selector from the model (including
    the parents).
    2. Slice the sequence.
    3. Update the model features locations using the crossmapper.

    TODO: Make sure everything is on the plus strand?

    :arg dict reference_model: Reference model.
    :arg str selector_id: Selector ID.
    :arg bool transcribe: Transcribe the sequence to RNA.
    :returns: RNA reference model.
    :rtype: dict
    """
    annotations = deepcopy(extract_feature_model(reference_model["annotations"], selector_id)[0])
    seq = (
        str(Seq(slice_to_selector(reference_model, selector_id)).transcribe()).lower()
        if transcribe
        else slice_to_selector(reference_model, selector_id)
    )
    selector_model = get_internal_selector_model(annotations, selector_id, True)

    rna_model = {"annotations": annotations, "sequence": {"seq": seq}}

    x = NonCoding(selector_model["exon"]).coordinate_to_noncoding
    g = Genomic().genomic_to_coordinate

    new_start = x(selector_model["exon"][0][0])[0] - 1
    new_end = x(selector_model["exon"][-1][-1])[0]

    for location, f_type in yield_locations(rna_model["annotations"]):
        if f_type in ["CDS", "exon"]:
            new_start_x = x(get_start(location))
            new_start_g = g(new_start_x[0])

            new_end_x = x(get_end(location))
            new_end_g = g(new_end_x[0] + new_end_x[1])

            set_start(location, new_start_g)
            set_end(location, new_end_g)
        else:
            set_start(location, new_start)
            set_end(location, new_end)
    return rna_model


def dna_to_rna(description):
    d = Description(description)
    d.to_delins()
    if d.errors:
        return {"errors": [splice_site([])]}

    delins = d.delins_model["variants"]
    sequences = d.get_sequences()
    exons = d.get_selector_model()["exon"]

    ref_seq = d.references["reference"]["sequence"]["seq"]
    alg_dna_variants, graph = extract_variants(ref_seq, delins_to_algebra(delins, sequences))
    local_supremal = get_local_supremal(ref_seq, graph)

    if (
        get_reference_mol_type(d.references["reference"]) == "genomic DNA"
        and _splice_sites_affected(exons, local_supremal)
    ):
        return {"errors": [splice_site([])]}

    alg_rna_sliced_variants = to_rna_variants(
        [algebra_variant_to_delins(v) for v in alg_dna_variants],
        d.get_sequences(),
        d.get_selector_model(),
    )
    rna_reference_models = get_rna_reference_models(d)
    rna_ref_seq = rna_reference_models["reference"]["sequence"]["seq"]
    alg_rna_variants, *_ = extract_variants(rna_ref_seq, delins_to_algebra(alg_rna_sliced_variants, {"reference": ref_seq})    )
    extracted_variants_model = to_model(to_hgvs(alg_rna_variants, rna_ref_seq), start_rule="variants")

    extracted_model = {
        "reference": d.corrected_model["reference"],
        "variants": extracted_variants_model,
        "predicted": True,
    }

    rna_selector_model = get_internal_selector_model(rna_reference_models["reference"]["annotations"], d.get_selector_id(), fix_exon=True)
    if rna_selector_model.get("cds"):
        x = Coding(rna_selector_model["exon"], rna_selector_model["cds"][0], d.is_inverted()).coordinate_to_coding
    else:
        x = NonCoding(rna_selector_model["exon"], d.is_inverted()).coordinate_to_noncoding
    for point, path in yield_point_locations_for_main_reference(extracted_model):
        set_by_path(extracted_model, path, coding_to_point(x(point["position"] - 1)))
    if d.is_inverted():
        for range_location, path in yield_ranges_main_reference(extracted_model):
            range_location["start"], range_location["end"] = range_location["end"], range_location["start"]
        for variant in extracted_model["variants"]:
            if variant.get("inserted"):
                for inserted in variant["inserted"]:
                    if inserted.get("sequence"):
                        inserted["sequence"] = reverse_complement(inserted["sequence"])
            if variant.get("deleted"):
                for deleted in variant["deleted"]:
                    if deleted.get("sequence"):
                        deleted["sequence"] = reverse_complement(deleted["sequence"])
        extracted_model["variants"].reverse()

    extracted_model["coordinate_system"] = "r"
    extracted_model["predicted"] = True

    return {"description": model_to_string(extracted_model)}


def rna_to_dna(description):
    d = Description(description)
    d.to_delins()
    if d.errors:
        return {"errors": [splice_site([])]}

    delins = d.delins_model["variants"]
    sequences = d.get_sequences()
    exons = d.get_selector_model()["exon"]

    ref_seq = d.references["reference"]["sequence"]["seq"]
    alg_dna_variants, graph = extract_variants(ref_seq, delins_to_algebra(delins, sequences))
    local_supremal = get_local_supremal(ref_seq, graph)

    if (
        get_reference_mol_type(d.references["reference"]) == "genomic DNA"
        and _splice_sites_affected(exons, local_supremal)
    ):
        return {"errors": [splice_site([])]}

    d.corrected_model["coordinate_system"] = get_coordinate_system_from_selector_id(
        d.references["reference"], d.get_selector_id()
    )

    return d.corrected_model
