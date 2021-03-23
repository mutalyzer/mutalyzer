from .description import Description
from .description_model import variant_to_description
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
    d.to_internal_indexing_model()
    sequences = d.get_sequences()
    output = []
    for i, variant in enumerate(d.internal_indexing_model["variants"]):
        if variant.get("location"):
            seq_del = slice_sequence(variant["location"], sequences["reference"])
        else:
            seq_del = ""
        seq_ins = get_inserted_sequence(variant, sequences)
        start = get_start(variant)
        end = get_end(variant)
        seq_l = sequences["reference"][start - left if start - left > 0 else 0 : start]
        seq_r = sequences["reference"][end : end + right]
        if variant.get("type") == "equal":
            output.append(
                {
                    "description": variant_to_description(
                        d.corrected_model["variants"][i]
                    ),
                    "left": seq_l,
                    "equal": shrink_seq(seq_del),
                    "right": seq_r,
                    "start": start - left if start - left > 0 else 0,
                }
            )
        else:
            output.append(
                {
                    "description": variant_to_description(
                        d.corrected_model["variants"][i]
                    ),
                    "left": seq_l,
                    "deleted": shrink_seq(seq_del),
                    "inserted": shrink_seq(seq_ins),
                    "right": seq_r,
                    "start": start - left if start - left > 0 else 0,
                }
            )
    return output
