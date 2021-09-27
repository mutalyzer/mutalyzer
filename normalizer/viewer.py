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


def view_variants(description, left=20, right=20):
    d = Description(description)
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
