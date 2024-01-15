from algebra import LCSgraph, Variant
from algebra.relations.graph_based import compare as compare_core

from mutalyzer import errors
from mutalyzer.description import Description
from mutalyzer.reference import retrieve_reference
from mutalyzer.util import get_end, get_inserted_sequence, get_start
from mutalyzer.viewer import view_delins


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
    reference_model = retrieve_reference(reference_id)[0]
    # TODO: update the reference id if different in the model?
    if reference_model is None:
        output["errors"] = [errors.reference_not_retrieved(reference_id, [])]
    else:
        output["sequence"] = reference_model["sequence"]["seq"]
        output["annotations"] = {"id": reference_model["annotations"]["id"]}
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
    if reference_type == "id":
        return _get_id(reference)


def _get_operand(m_input, m_type, reference):
    if m_type == "sequence":
        return _get_sequence(m_input)
    if m_type == "hgvs":
        return _get_hgvs_and_variant(m_input)
    if m_type == "variant":
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


def _algebra_variant_to_delins(variant):
    delins_variant = {
        "type": "deletion_insertion",
        "source": "reference",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": variant.start},
            "end": {"type": "point", "position": variant.end},
        },
        "deleted": [],
        "inserted": [],
    }
    if variant.sequence:
        delins_variant["inserted"].append(
            {"sequence": variant.sequence, "source": "description"}
        )

    return delins_variant


def _algebra_variant_to_name_model(variant):
    def _position_to_hgvs():
        if variant.end - variant.start == 1:
            return {"type": "point", "position": variant.start + 1}
        if variant.start == variant.end:
            return {
                "type": "range",
                "start": {"type": "point", "position": variant.start},
                "end": {"type": "point", "position": variant.start + 1},
            }
        return {
            "type": "range",
            "start": {"type": "point", "position": variant.start + 1},
            "end": {"type": "point", "position": variant.end},
        }

    delins_variant = {
        "source": "reference",
        "location": _position_to_hgvs(),
        "deleted": [],
        "inserted": [],
    }
    if variant.sequence:
        if variant.start == variant.end:
            delins_variant["type"] = "insertion"
        else:
            delins_variant["type"] = "deletion_insertion"
        delins_variant["inserted"].append(
            {"sequence": variant.sequence, "source": "description"}
        )
    else:
        delins_variant["type"] = "deletion"
    return delins_variant


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

    lhs_alg_variants = _get_algebra_variants(lhs_d)
    lhs_graph = LCSgraph.from_variant(lhs_reference, lhs_alg_variants)
    lhs_supremal = lhs_graph.supremal

    rhs_alg_variants = _get_algebra_variants(rhs_d)
    rhs_graph = LCSgraph.from_variant(rhs_reference, rhs_alg_variants)
    rhs_supremal = rhs_graph.supremal

    output["relation"] = compare_core(lhs_reference, lhs_graph, rhs_graph).value

    if lhs_d.corrected_model.get("reference"):
        ref_id = lhs_d.corrected_model["reference"]["id"]
    elif lhs_d.references["reference"].get("annotations") and lhs_d.references[
        "reference"
    ]["annotations"].get("id"):
        ref_id = lhs_d.references["reference"]["annotations"]["id"]
    else:
        ref_id = rhs_reference

    output["supremal_lhs"] = {
        "hgvs": f"{ref_id}:g.{lhs_supremal.to_hgvs()}",
        "spdi": lhs_supremal.to_spdi(ref_id),
    }

    output["supremal_rhs"] = {
        "hgvs": f"{ref_id}:g.{rhs_supremal.to_hgvs()}",
        "spdi": rhs_supremal.to_spdi(ref_id),
    }

    lhs_supremal_delins = [_algebra_variant_to_delins(lhs_supremal)]
    rhs_supremal_delins = [_algebra_variant_to_delins(rhs_supremal)]
    output["view_lhs_supremal"] = view_delins(
        lhs_supremal_delins,
        [_algebra_variant_to_name_model(lhs_supremal)],
        lhs_d.get_sequences(),
    )
    output["view_rhs_supremal"] = view_delins(
        rhs_supremal_delins,
        [_algebra_variant_to_name_model(rhs_supremal)],
        rhs_d.get_sequences(),
    )

    return output


def compare_sequences_based(reference, reference_type, lhs, lhs_type, rhs, rhs_type):
    output = {}

    c_reference = _get_reference(reference, reference_type)
    if c_reference.get("errors"):
        _extend_errors(output, "reference", c_reference["errors"])
        return output
    ref_seq = c_reference["sequence"]
    c_lhs = _get_operand(lhs, lhs_type, ref_seq)
    c_rhs = _get_operand(rhs, rhs_type, ref_seq)
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

    lhs_graph = LCSgraph.from_sequence(ref_seq, lhs_seq)
    lhs_supremal = lhs_graph.supremal

    rhs_graph = LCSgraph.from_sequence(ref_seq, rhs_seq)
    rhs_supremal = rhs_graph.supremal

    output["relation"] = compare_core(ref_seq, lhs_graph, rhs_graph).value

    output["supremal_lhs"] = {
        "hgvs": f"{lhs_supremal.to_hgvs()}",
        "spdi": lhs_supremal.to_spdi(),
    }

    output["supremal_rhs"] = {
        "hgvs": f"{rhs_supremal.to_hgvs()}",
        "spdi": rhs_supremal.to_spdi(),
    }

    lhs_supremal_delins = [_algebra_variant_to_delins(lhs_supremal)]
    rhs_supremal_delins = [_algebra_variant_to_delins(rhs_supremal)]
    output["view_lhs_supremal"] = view_delins(
        lhs_supremal_delins,
        [_algebra_variant_to_name_model(lhs_supremal)],
        {"reference": ref_seq},
    )
    output["view_rhs_supremal"] = view_delins(
        rhs_supremal_delins,
        [_algebra_variant_to_name_model(rhs_supremal)],
        {"reference": ref_seq},
    )

    return output


def compare_hgvs_based(reference, reference_type, lhs, lhs_type, rhs, rhs_type):
    """
    Compare two HGVS descriptions (either complete, i.e., including the
    reference id, or just the variants relative to reference sequence indicated
    by an id or directly as a string).
    """
    if lhs_type == "hgvs" and rhs_type == "hgvs":
        lhs_d = Description(lhs)
        lhs_d.to_delins()
        rhs_d = Description(rhs)
        rhs_d.to_delins()
    elif lhs_type == "hgvs" and rhs_type == "variant":
        lhs_d = Description(lhs)
        lhs_d.to_delins()
        if lhs_d.get_sequences() and lhs_d.get_sequences().get("reference"):
            rhs_d = Description(
                description=rhs,
                only_variants=True,
                sequence=lhs_d.get_sequences()["reference"],
            )
            rhs_d.to_delins()
        else:
            return
    elif lhs_type == "variant" and rhs_type == "variant":
        if reference_type == "sequence":
            reference_sequence = reference
            ref_id = reference
        elif reference_type == "id":
            check = _get_id(reference)
            if check.get("errors"):
                return check
            reference_sequence = check["sequence"]
            ref_id = check["annotations"]["id"]
        lhs_d = Description(
            description=lhs, only_variants=True, sequence=reference_sequence
        )
        lhs_d.to_delins()
        if not lhs_d.errors:
            lhs_d.references["reference"]["annotations"] = {"id": ref_id}

        rhs_d = Description(
            description=rhs, only_variants=True, sequence=reference_sequence
        )
        rhs_d.to_delins()
        if not rhs_d.errors:
            rhs_d.references["reference"]["annotations"] = {"id": ref_id}
    return compare_hgvs(lhs_d, rhs_d)


def compare(reference, reference_type, lhs, lhs_type, rhs, rhs_type):
    """
    Generic interface to the algebra.

    Parameters
    ----------
    reference : str
        The reference id or the sequence.
    reference_type : str
        If the reference is a sequence or an id.
    lhs : str
        The left operand.
    lhs_type : str
        The type of the left operand.
    rhs : str
        The right operand.
    rhs_type : str
        The type of the right operand.

    Returns
    -------
    dict
        The relation and additional information.
    """
    checks = _input_types_check(reference_type, lhs_type, rhs_type)

    if checks:
        return checks

    if reference_type in ["sequence", "id"] and (
        lhs_type == "sequence" or rhs_type == "sequence"
    ):
        return compare_sequences_based(
            reference, reference_type, lhs, lhs_type, rhs, rhs_type
        )
    return compare_hgvs_based(reference, reference_type, lhs, lhs_type, rhs, rhs_type)
