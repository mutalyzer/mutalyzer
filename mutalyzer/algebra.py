from algebra import Variant
from algebra.lcs import edit, lcs_graph
from algebra.relations.sequence_based import compare as compare_core
from algebra.relations.supremal_based import compare as compare_supremal
from algebra.relations.supremal_based import find_supremal, spanning_variant
from algebra.variants import Variant

import mutalyzer.errors as errors
from mutalyzer.description import Description
from mutalyzer.reference import retrieve_reference
from mutalyzer.util import get_end, get_inserted_sequence, get_start
from mutalyzer.viewer import view_variants, view_variants_normalized


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
    reference_model = retrieve_reference(reference_id)[0]
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


def _influence_interval_supremal(supremal):
    min_pos = supremal.start
    max_pos = supremal.end
    if supremal == Variant(0, 0, ""):
        return {"equal": True}
    else:
        return {"min_pos": min_pos, "max_pos": max_pos}


def _influence_interval(output, ref_seq, obs_seq, hs):
    _, lcs_nodes = edit(ref_seq, obs_seq)
    _, edges = lcs_graph(ref_seq, obs_seq, lcs_nodes)
    output[f"influence_{hs}"] = _influence_interval_supremal(
        spanning_variant(ref_seq, obs_seq, edges)
    )


def _check_sequences_equality(output, lhs, rhs):
    if lhs["reference_sequence"] != rhs["reference_sequence"]:
        _append_error(
            output,
            "reference",
            {
                "code": "ESEQUENCEMISMATCH",
                "details": "Different reference sequences for LHS and RHS.",
            },
        )


def _check_sequence_length(output, ref_seq, len_max):
    print(len(ref_seq))
    if len(ref_seq) > len_max:
        _append_error(
            output,
            "reference",
            {
                "code": "ESEQUENCELENGTH",
                "details": f"Sequence length {len(ref_seq)} too large (maximum supported is {len_max}).",
            },
        )


def _get_algebra_variants(description):
    """
    Convert the delins variants to algebra variant edges.

    :param description: Normalizer description object.
    :return: Algebra variant edges.
    """
    edges = []
    for variant in description.delins_model["variants"]:
        edges.append(
            Variant(
                get_start(variant),
                get_end(variant),
                get_inserted_sequence(variant, description.get_sequences()),
            )
        )
    return edges


def compare_hgvs(lhs_d, rhs_d):
    output = {}

    if lhs_d.errors:
        _extend_errors(output, "lhs", lhs_d.errors)

    if rhs_d.errors:
        _extend_errors(output, "rhs", rhs_d.errors)

    if output:
        return output

    lhs_reference = lhs_d.get_sequences()["reference"]
    rhs_reference = rhs_d.get_sequences()["reference"]

    if lhs_reference != rhs_reference:
        _append_error(
            output,
            "reference",
            {
                "code": "ESEQUENCEMISMATCH",
                "details": "Different reference sequences for LHS and RHS.",
            },
        )

    if output:
        return output

    lhs_observed = lhs_d.get_sequences()["observed"]
    lhs_algebra_variants = _get_algebra_variants(lhs_d)
    lhs_spanning = spanning_variant(lhs_reference, lhs_observed, lhs_algebra_variants)
    lhs_supremal = find_supremal(lhs_reference, lhs_spanning)

    rhs_observed = rhs_d.get_sequences()["observed"]
    rhs_algebra_variants = _get_algebra_variants(rhs_d)
    rhs_spanning = spanning_variant(rhs_reference, rhs_observed, rhs_algebra_variants)
    rhs_supremal = find_supremal(rhs_reference, rhs_spanning)

    output["relation"] = compare_supremal(
        lhs_reference, lhs_supremal, rhs_supremal
    ).value

    output["influence_lhs"] = _influence_interval_supremal(lhs_supremal)
    output["influence_rhs"] = _influence_interval_supremal(rhs_supremal)

    output["view_lhs"] = view_variants_normalized(lhs_d)
    output["view_rhs"] = view_variants_normalized(rhs_d)

    return output


def compare_sequences_based(reference, reference_type, lhs, lhs_type, rhs, rhs_type):
    output = {}

    c_reference = _get_reference(reference, reference_type)
    if c_reference.get("errors"):
        _extend_errors(output, "reference", c_reference["errors"])
        return output
    ref_seq = c_reference["sequence"]
    c_lhs = _get_operator(lhs, lhs_type, ref_seq)
    c_rhs = _get_operator(rhs, rhs_type, ref_seq)
    c_lhs["reference_sequence"] = ref_seq
    c_rhs["reference_sequence"] = ref_seq

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

    output["relation"] = compare_core(ref_seq, lhs_seq, rhs_seq).value
    _influence_interval(output, ref_seq, lhs_seq, "lhs")
    _influence_interval(output, ref_seq, rhs_seq, "rhs")

    if c_lhs.get("view"):
        output["view_lhs"] = c_lhs["view"]
    if c_rhs.get("view"):
        output["view_rhs"] = c_rhs["view"]

    return output


def compare_hgvs_based(reference, reference_type, lhs, lhs_type, rhs, rhs_type):
    if lhs_type == "hgvs" and rhs_type == "hgvs":
        lhs_d = Description(lhs)
        lhs_d.normalize()
        rhs_d = Description(rhs)
        rhs_d.normalize()
    elif lhs_type == "hgvs" and rhs_type == "variant":
        lhs_d = Description(lhs)
        lhs_d.normalize()
        if lhs_d.get_sequences() and lhs_d.get_sequences().get("reference"):
            rhs_d = Description(
                description=rhs,
                only_variants=True,
                sequence=lhs_d.get_sequences()["reference"],
            )
            rhs_d.normalize()
        else:
            return
    elif lhs_type == "variant" and rhs_type == "variant":
        if reference_type == "sequence":
            reference_sequence = "reference"
        elif reference_type == "id":
            check = _get_id(reference)
            if check.get("errors"):
                return check
            else:
                reference_sequence = check["sequence"]
        lhs_d = Description(
            description=lhs, only_variants=True, sequence=reference_sequence
        )
        lhs_d.normalize()
        rhs_d = Description(
            description=rhs, only_variants=True, sequence=reference_sequence
        )
        rhs_d.normalize()
    return compare_hgvs(lhs_d, rhs_d)


def compare(reference, reference_type, lhs, lhs_type, rhs, rhs_type):
    checks = _input_types_check(reference_type, lhs_type, rhs_type)

    if checks:
        return checks

    if reference_type in ["sequence", "id"] and (
        lhs_type == "sequence" or rhs_type == "sequence"
    ):
        return compare_sequences_based(
            reference, reference_type, lhs, lhs_type, rhs, rhs_type
        )
    else:
        return compare_hgvs_based(
            reference, reference_type, lhs, lhs_type, rhs, rhs_type
        )
