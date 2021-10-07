import copy

from mutalyzer_mutator.mutator import reverse_complement

from .converter.to_delins import variant_to_delins
from .description import Description
from .description_model import variant_to_description, variants_to_description
from .util import get_end, get_inserted_sequence, get_start, slice_sequence


def shrink_seq(seq, limit=100, left=50, right=50):
    if len(seq) > limit:
        return {
            "left": seq[:left],
            "right": seq[-right:],
            "original_length": len(seq),
            "shrunk": True,
        }
    else:
        return {"seq": seq}


def view_variants_old(
    description, only_variants=False, sequence=None, left=20, right=20
):
    d = Description(
        description=description, only_variants=only_variants, sequence=sequence
    )
    d.normalize()
    if d.errors:
        return d.output()
    sequences = d.get_sequences()
    output = []
    for i, variant in enumerate(d.internal_indexing_model["variants"]):
        delins_variant = variant_to_delins(variant)
        if delins_variant is None:
            delins_variant = copy.deepcopy(variant)
        if delins_variant.get("location"):
            seq_del = slice_sequence(delins_variant["location"], sequences["reference"])
        else:
            seq_del = ""
        seq_ins = get_inserted_sequence(delins_variant, sequences)

        v_start = get_start(delins_variant)
        v_end = get_end(delins_variant)
        start = v_start - left if v_start - left > 0 else 0
        seq_left = sequences["reference"][start:v_start]
        seq_right = sequences["reference"][v_end : v_end + right]
        details = {
            "description": variant_to_description(d.corrected_model["variants"][i]),
            "seq_length": len(sequences["reference"]),
        }
        if d.is_inverted():
            seq_left, seq_right = reverse_complement(seq_right), reverse_complement(
                seq_left
            )
            seq_del = reverse_complement(seq_del)
            seq_ins = reverse_complement(seq_ins)
            details["description"] = variant_to_description(
                d.corrected_model["variants"][-1 - i]
            )
            details["inverted"] = True
            start = (
                v_end - 1 + right
                if v_end + right < len(sequences["reference"])
                else len(sequences["reference"]) - 1
            )

        if variant.get("type") in [None, "equal"]:
            details.update(
                {
                    "left": seq_left,
                    "equal": shrink_seq(seq_del),
                    "right": seq_right,
                    "start": start,
                }
            )
        else:
            details.update(
                {
                    "left": seq_left,
                    "deleted": shrink_seq(seq_del),
                    "inserted": shrink_seq(seq_ins),
                    "right": seq_right,
                    "start": start,
                }
            )
        output.append(details)
    if d.is_inverted():
        output.reverse()
    return output


def _get_sequence_view(seq, l_l=10, l_r=10):
    s = 0
    e = len(seq)
    if e - s > l_l + l_r:
        return {"left": seq[s : s + l_l], "right": seq[e - l_l : e]}
    if 0 < e - s <= l_l + l_r:
        return {"sequence": seq[s:e]}


def _get_view_inside(s, e, ref_seq, variant):
    view = {"start": s, "end": e, "type": "variant"}
    del_seq = ref_seq[s:e]
    if del_seq:
        view["deleted"] = _get_sequence_view(del_seq)
    ins_seq = get_inserted_sequence(variant, {"reference": ref_seq})
    if ins_seq:
        view["inserted"] = _get_sequence_view(ins_seq)
        view["inserted"]["length"] = len(ins_seq)
    return view


def _get_view_outside(s, e, ref_seq):
    view = {"start": s, "end": e, "type": "outside"}
    seq = ref_seq[s:e]
    if seq:
        view.update(_get_sequence_view(seq))
    return view


def _get_segments(variants, ref_seq):
    points = []
    for variant in variants:
        points += [get_start(variant), get_end(variant)]
    points = [0] + points + [len(ref_seq)]
    return [(points[i - 1], points[i]) for i in range(len(points))[1:]]


def _invert_views(views, ref_length):
    inv_views = []
    for view in views[::-1]:
        inv_view = copy.deepcopy(view)
        inv_view["start"] = ref_length - view["end"]
        inv_view["end"] = ref_length - view["start"]
        if inv_view.get("type") == "outside":
            if inv_view.get("left"):
                inv_view["left"] = reverse_complement(view["right"])
                inv_view["right"] = reverse_complement(view["left"])
            elif inv_view.get("sequence"):
                inv_view["sequence"] = reverse_complement(view["sequence"])
        if inv_view.get("type") == "variant":
            if inv_view.get("deleted") and inv_view["deleted"].get("left"):
                inv_view["deleted"]["left"] = reverse_complement(
                    view["deleted"]["right"]
                )
                inv_view["deleted"]["right"] = reverse_complement(
                    view["deleted"]["left"]
                )
            elif inv_view.get("deleted") and inv_view["deleted"].get("sequence"):
                inv_view["deleted"]["sequence"] = reverse_complement(view["deleted"]["sequence"])
        inv_views.append(inv_view)
    return inv_views


def view_variants(
    description, only_variants=False, sequence=None, left=10, right=10, invert=True
):
    d = Description(
        description=description, only_variants=only_variants, sequence=sequence
    )
    d.normalize()
    if d.errors:
        return d.output()

    ref_seq = d.get_sequences()["reference"]

    segments = _get_segments(d.delins_model["variants"], ref_seq)
    views = []
    for i, segment in enumerate(segments):
        if i % 2 == 0:
            view = _get_view_outside(*segment, ref_seq)
        else:
            view = {
                "description": variant_to_description(
                    d.corrected_model["variants"][i // 2]
                )
            }
            view.update(
                _get_view_inside(*segment, ref_seq, d.delins_model["variants"][i // 2])
            )
        views.append(view)
    if invert and d.is_inverted():
        return {
            "views": _invert_views(views, len(ref_seq)),
            "seq_length": len(ref_seq),
            "inverted": True,
        }
    return {"views": views, "seq_length": len(ref_seq)}
