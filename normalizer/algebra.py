from algebra.algebra import compare as compare_core
from mutalyzer_hgvs_parser import to_model
from mutalyzer_hgvs_parser.exceptions import UnexpectedCharacter, UnexpectedEnd
from mutalyzer_mutator import mutate as mutalyzer_mutator

import normalizer.errors as errors
from normalizer.description import Description
from normalizer.reference import retrieve_reference

from .converter.to_delins import to_delins
from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing


def _get_variant(variant, reference):
    output = {"input": variant, "type": "variant"}

    try:
        variants_model = to_model(variant, "variants")
    except UnexpectedCharacter as e:
        output["errors"] = [errors.syntax_uc(e)]
        return output
    except UnexpectedEnd as e:
        output["errors"] = [errors.syntax_ueof(e)]
        return output

    artificial_model = {
        "coordinate_system": "g",
        "variants": variants_model,
    }
    sequence = {"sequence": {"seq": reference["sequence"]}}
    delins_model = to_delins(
        to_internal_indexing(to_internal_coordinates(artificial_model, sequence))
    )
    output["sequence"] = mutalyzer_mutator(
        {"reference": reference["sequence"]}, delins_model["variants"]
    )
    return output


def _get_hgvs(description):
    output = {"input": description, "type": "hgvs"}
    d = Description(description)

    d.normalize()
    status = d.output()
    if status.get("input_model"):
        output["reference"] = status["input_model"]["reference"]
    if status.get("errors"):
        output["errors"] = status["errors"]
        return output

    if status.get("infos"):
        output["infos"] = status["infos"]

    sequences = d.get_sequences()
    output["sequence"] = {
        "reference": sequences["reference"],
        "observed": sequences["observed"],
    }

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
    else:
        return {
            "errors": [
                {
                    "code": "ENOTSUPPORTED",
                    "details": f"'{reference_type}' not supported as input type.",
                }
            ]
        }


def _get_operator(m_input, m_type, reference):
    if m_type == "sequence":
        return _get_sequence(m_input)
    elif m_type == "hgvs":
        return _get_hgvs(m_input)
    elif m_type == "variant":
        return _get_variant(m_input, reference)
    else:
        return {
            "errors": [
                {
                    "code": "ENOTSUPPORTED",
                    "details": f"'{m_type}' not supported as input type.",
                }
            ]
        }


def compare(reference, reference_type, lhs, lhs_type, rhs, rhs_type):
    ref_seq = ""

    c_reference = _get_reference(reference, reference_type)
    c_lhs = _get_operator(lhs, lhs_type, c_reference)
    c_rhs = _get_operator(rhs, rhs_type, c_reference)

    if c_reference and c_reference.get("errors"):
        return {"errors": c_reference["errors"]}
    if c_reference:
        ref_seq = c_reference["sequence"]
    obs_seq_1 = c_lhs["sequence"]
    obs_seq_2 = c_rhs["sequence"]

    return {"relation": compare_core(ref_seq, obs_seq_1, obs_seq_2)[0]}
