from copy import deepcopy

from mutalyzer_mutator.mutator import reverse_complement

from .description import Description
from .description_model import variant_to_description
from .util import get_end, get_inserted_sequence, get_start


def _get_sequence_view(seq, l_l=15, l_r=15):
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


def _invert_left_right(view, inv_view):
    inv_view["left"] = reverse_complement(view["right"])
    inv_view["right"] = reverse_complement(view["left"])


def _invert_view(v, inv_v):
    if inv_v.get("left") and inv_v.get("right"):
        _invert_left_right(v, inv_v)
    elif inv_v.get("sequence"):
        inv_v["sequence"] = reverse_complement(v["sequence"])


def _invert_views(views, ref_length):
    inv_vs = []
    for v in views[::-1]:
        inv_v = deepcopy(v)
        inv_v["start"] = ref_length - v["end"]
        inv_v["end"] = ref_length - v["start"]
        if inv_v.get("type") == "outside":
            _invert_view(v, inv_v)
        if inv_v.get("type") == "variant":
            if inv_v.get("deleted"):
                _invert_view(v["deleted"], inv_v["deleted"])
            if inv_v.get("inserted"):
                _invert_view(v["inserted"], inv_v["inserted"])
        inv_vs.append(inv_v)
    return inv_vs


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
