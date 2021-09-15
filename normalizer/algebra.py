from algebra.algebra import compare as compare_core

import normalizer.errors as e
from normalizer.reference import retrieve_reference

from .mutator import mutate, mutate_sequence


def compare(reference, reference_type, lhs, lhs_type, rhs, rhs_type):
    ref_seq = ""
    obs_seq_1 = ""
    obs_seq_2 = ""
    errors = []

    if reference_type == "sequence":
        # TODO: Check if it is a DNA sequence?
        ref_seq = reference
    elif reference_type == "id":
        reference_model = retrieve_reference(reference)
        if reference_model is None:
            errors.append(e.reference_not_retrieved(reference, ["reference"]))
        else:
            ref_seq = reference_model["sequence"]["seq"]

    if lhs_type == "sequence":
        obs_seq_1 = lhs
    elif lhs_type == "hgvs":
        mutated = mutate(lhs)
        if mutated.get("errors"):
            errors.extend(mutated["errors"])
        else:
            obs_seq_1 = mutated["sequence"]["seq"]
    elif lhs_type == "variant":
        obs_seq_1 = mutate_sequence(ref_seq, lhs)

    if rhs_type == "sequence":
        obs_seq_2 = rhs
    elif rhs_type == "hgvs":
        mutated = mutate(rhs)
        if mutated.get("errors"):
            errors.extend(mutated["errors"])
        else:
            obs_seq_2 = mutated["sequence"]["seq"]
    elif lhs_type == "variant":
        obs_seq_2 = mutate_sequence(ref_seq, rhs)
    if not errors:
        return {"relation": compare_core(ref_seq, obs_seq_1, obs_seq_2)[0]}
    else:
        return {"errors": errors}
