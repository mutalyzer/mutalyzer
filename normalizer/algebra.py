from algebra.algebra import compare as compare_core
from algebra.influence_interval import influence_interval

import normalizer.errors as errors
from normalizer.description import Description
from normalizer.reference import retrieve_reference
from normalizer.description_extractor import description_extractor
from .util import get_end, get_inserted_sequence, get_start, slice_sequence
from normalizer.description_model import variant_to_description


def _get_hgvs_and_variant(variant, only_variants=False, ref_seq=None):
    output = {"input": variant, "type": "variant" if only_variants else "hgvs"}

    d = Description(description=variant, only_variants=only_variants, sequence=ref_seq)
    d.normalize()

    status = d.output()
    if status.get("input_model") and status["input_model"].get("reference"):
        output["reference"] = status["input_model"]["reference"]
    if status.get("errors"):
        output["errors"] = status["errors"]
        return output
    if status.get("infos"):
        output["infos"] = status["infos"]

    sequences = d.get_sequences()
    output["sequence"] = sequences["observed"]
    output["reference_sequence"] = sequences["reference"]

    return output


def _get_id(reference_id):
    output = {"input": reference_id, "type": "id", "reference": {"id": reference_id}}
    reference_model = retrieve_reference(reference_id)
    # TODO: update the reference id if different in the model?
    if reference_model is None:
        output["errors"] = [errors.reference_not_retrieved(reference_id, [])]
    else:
        output["sequence"] = reference_model["sequence"]["seq"]
    return output


def _get_sequence(seq):
    return {"input": seq, "type": "sequence", "sequence": seq}


def _get_reference(reference, reference_type):
    if reference is None and reference_type is None:
        return None
    if reference_type == "sequence":
        return _get_sequence(reference)
    elif reference_type == "id":
        return _get_id(reference)


def _get_operator(m_input, m_type, reference):
    if m_type == "sequence":
        return _get_sequence(m_input)
    elif m_type == "hgvs":
        return _get_hgvs_and_variant(m_input)
    elif m_type == "variant":
        return _get_hgvs_and_variant(m_input, True, reference)


def _sort(d):
    if isinstance(d, list):
        for v in d:
            _sort(v)
    if isinstance(d, dict):
        for k in d:
            if isinstance(d[k], list):
                d[k] = sorted(d[k])
            if isinstance(d[k], dict):
                _sort(d[k])


def _add_message(o, m, t):
    if m and m.get("errors"):
        _sort(m["errors"])
        if o.get("errors") is None:
            o["errors"] = {}
        o["errors"][t] = m["errors"]


def _add_standard_messages(o, c_reference, c_lhs, c_rhs):
    _add_message(o, c_reference, "reference")
    _add_message(o, c_lhs, "lhs")
    _add_message(o, c_rhs, "rhs")


def _append_error(d, k, v):
    if d.get("errors") is None:
        d["errors"] = {}
    if d["errors"].get(k) is None:
        d["errors"][k] = []
    d["errors"][k].append(v)


def _extend_errors(d, k, v):
    if d.get("errors") is None:
        d["errors"] = {}
    if d["errors"].get(k) is None:
        d["errors"][k] = []
    d["errors"][k].extend(v)


def _input_types_check(reference_type, lhs_type, rhs_type):
    output = {}
    if reference_type and reference_type not in ["sequence", "id"]:
        _append_error(
            output,
            "reference_type",
            errors.invalid_input(reference_type, ["sequence", "id"]),
        )

    operator_types = ["sequence", "variant", "hgvs"]
    if lhs_type not in operator_types:
        _append_error(
            output,
            "lhs_type",
            errors.invalid_input(lhs_type, operator_types),
        )

    if rhs_type not in operator_types:
        _append_error(
            output,
            "rhs_type",
            errors.invalid_input(rhs_type, operator_types),
        )
    return output


def _influence_interval(output, ref_seq, obs_seq, hs):
    try:
        min_pos, max_pos = influence_interval(ref_seq, obs_seq)
    except Exception as e:
        if "Equal to the reference" in str(e):
            output[f"influence_{hs}"] = {"equal": True}
    else:
        output[f"influence_{hs}"] = {"min_pos": min_pos, "max_pos": max_pos}


def _add_sequences(output, ref_seq, lhs_seq, rhs_seq, len_max=100):
    if len(ref_seq) <= len_max:
        output["ref_seq"] = ref_seq
    if len(lhs_seq) <= len_max:
        output["lhs_seq"] = lhs_seq
    if len(rhs_seq) <= len_max:
        output["rhs_seq"] = rhs_seq


def _influence(ref_seq, obs_seq):
    influence = {}
    try:
        min_pos, max_pos = influence_interval(ref_seq, obs_seq)
    except Exception as e:
        if "Equal to the reference" in str(e):
            influence = {"equal": True}
    else:
        influence = {"min_pos": min_pos, "max_pos": max_pos}


def _get_segments(variants, ref_seq):
    points = []
    for variant in variants:
        points += [get_start(variant), get_end(variant)]
    points = [0] + points + [len(ref_seq)]
    return [(points[i - 1], points[i]) for i in range(len(points))[1:]]


def _get_sequence_view(seq, l_l=2, l_r=2):
    s = 0
    e = len(seq)
    if e - s > l_l + l_r:
        return {"left": seq[s : s + l_l], "right": seq[e - l_l : e]}
    if 0 < e - s <= l_l + l_r:
        return {"sequence": seq[s:e]}


def _get_view_outside(s, e, ref_seq, l_l=2, l_r=2):
    view = {"start": s, "end": e, "type": "outside"}
    if e - s > l_l + l_r:
        view["left"] = ref_seq[s : s + l_l]
        view["right"] = ref_seq[e - l_l : e]
    if 0 < e - s <= l_l + l_r:
        view["sequence"] = ref_seq[s:e]
    return view


def _get_view_inside(s, e, ref_seq, variant, l_l=2, l_r=2):
    view = {"start": s, "end": e, "type": "variant"}
    if e - s > l_l + l_r:
        view["deleted"] = {"left": ref_seq[s : s + l_l], "right": ref_seq[e - l_l : e]}
    if 0 < e - s <= l_l + l_r:
        view["deleted"] = {"sequence": ref_seq[s:e]}
    ins_seq = get_inserted_sequence(variant, {"reference": ref_seq})
    if ins_seq:
        view["inserted"] = _get_sequence_view(ins_seq)
    return view


def details(ref_seq, obs_seq):
    # influence = _influence(ref_seq, obs_seq)
    variants = description_extractor(ref_seq, obs_seq)
    d = Description(variants, only_variants=True, sequence=ref_seq)
    d.normalize()
    segments = _get_segments(d.delins_model["variants"], ref_seq)
    view = []
    for i, segment in enumerate(segments):
        if i % 2 == 0:
            view.append(_get_view_outside(*segment, ref_seq))
        else:
            view.append(
                _get_view_inside(*segment, ref_seq, d.delins_model["variants"][i // 2])
            )
    return view


def add_details():
    pass


def compare(reference, reference_type, lhs, lhs_type, rhs, rhs_type):
    output = _input_types_check(reference_type, lhs_type, rhs_type)
    if output:
        return output

    c_reference = _get_reference(reference, reference_type)

    if reference_type and reference is None:
        _append_error(output, "reference", errors.missing_parameter("reference"))
        return output

    if reference_type in ["sequence", "id"]:
        c_reference = _get_reference(reference, reference_type)
        if c_reference.get("errors"):
            _extend_errors(output, "reference", c_reference["errors"])
            return output
        ref_seq = c_reference["sequence"]
        c_lhs = _get_operator(lhs, lhs_type, ref_seq)
        c_rhs = _get_operator(rhs, rhs_type, ref_seq)
    elif lhs_type == "hgvs":
        c_lhs = _get_hgvs_and_variant(lhs)
        if c_lhs.get("errors"):
            _extend_errors(output, "lhs", c_lhs["errors"])
            return output
        c_rhs = _get_operator(rhs, rhs_type, c_lhs["reference_sequence"])

    _add_standard_messages(output, c_reference, c_lhs, c_rhs)

    if output.get("errors"):
        return output

    if c_reference:
        ref_seq = c_reference["sequence"]
    else:
        ref_seq = c_lhs["reference_sequence"]
    lhs_seq = c_lhs["sequence"]
    rhs_seq = c_rhs["sequence"]

    len_max = 100000
    if len(ref_seq) > len_max:
        _append_error(
            output,
            "reference",
            {
                "code": "ESEQUENCELENGTH",
                "details": f"Sequence length {len(ref_seq)} too large (maximum supported is {len_max}).",
            },
        )

    if output.get("errors"):
        return output

    output["relation"] = compare_core(ref_seq, lhs_seq, rhs_seq)[0]
    _influence_interval(output, ref_seq, lhs_seq, "lhs")
    _influence_interval(output, ref_seq, rhs_seq, "rhs")
    _add_sequences(output, ref_seq, lhs_seq, rhs_seq)

    output["view_lhs"] = details(ref_seq, lhs_seq)
    output["view_rhs"] = details(ref_seq, rhs_seq)

    return output
