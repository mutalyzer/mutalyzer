from algebra.algebra import compare as compare_core
from algebra.influence_interval import influence_interval

import mutalyzer.errors as errors
from mutalyzer.description import Description
from mutalyzer.reference import retrieve_reference
from mutalyzer.viewer import view_variants


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
    output["view"] = view_variants(
        description=variant, only_variants=only_variants, sequence=ref_seq
    )

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
    return {
        "input": seq,
        "type": "sequence",
        "sequence": seq,
        "reference_sequence": seq,
    }


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


def _check_sequences_equality(output, lhs, rhs):
    if lhs["reference_sequence"] != rhs["reference_sequence"]:
        _append_error(
            output,
            "reference",
            {
                "code": "ESEQUENCEMISMATCH",
                "details": f"Different reference sequences for LHS and RHS.",
            },
        )


def _check_sequence_length(output, ref_seq, len_max):
    if len(ref_seq) > len_max:
        _append_error(
            output,
            "reference",
            {
                "code": "ESEQUENCELENGTH",
                "details": f"Sequence length {len(ref_seq)} too large (maximum supported is {len_max}).",
            },
        )


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
        c_lhs["reference_sequence"] = ref_seq
        c_rhs["reference_sequence"] = ref_seq
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

    _check_sequence_length(output, ref_seq, 100000)
    _check_sequences_equality(output, c_lhs, c_rhs)

    if output.get("errors"):
        return output

    output["relation"] = compare_core(ref_seq, lhs_seq, rhs_seq)[0]
    _influence_interval(output, ref_seq, lhs_seq, "lhs")
    _influence_interval(output, ref_seq, rhs_seq, "rhs")

    if c_lhs.get("view"):
        output["view_lhs"] = c_lhs["view"]
    if c_rhs.get("view"):
        output["view_rhs"] = c_rhs["view"]

    return output
