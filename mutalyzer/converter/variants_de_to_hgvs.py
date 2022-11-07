import copy

from mutalyzer.util import (
    get_end,
    get_inserted_sequence,
    get_location_length,
    get_start,
    roll,
    slice_sequence,
)


def seq2repeats(long_sequence):
    for i in range(1, len(long_sequence) + 1):
        for j in range(1, len(long_sequence) + 1):
            if i * j == len(long_sequence):
                if long_sequence[:i] * j == long_sequence:
                    yield long_sequence[:i], j
                break


def seq_present_before(observed, ins_seq, start, end):
    for repeat, i in seq2repeats(ins_seq):
        if observed[start - len(repeat) : end] == repeat:
            return repeat, i
    return "", 0


def is_deletion(delins_variant):
    if (delins_variant.get("inserted") is None) or (
        len(delins_variant.get("inserted")) == 0
    ):
        return True
    for inserted in delins_variant["inserted"]:
        if get_location_length(inserted["location"]) > 0:
            return False
    return True


def update_inserted_with_sequences(inserted, sequences):
    new_inserted = []
    for insert in inserted:
        if insert["source"] == "observed":
            seq = sequences["observed"][
                get_start(insert["location"]) : get_end(insert["location"])
            ]
            if seq:
                new_inserted.append({"sequence": seq, "source": "description"})
        else:
            new_inserted.append(insert)
    return new_inserted


def de_variants_clean(variants, sequences=None):
    """
    Apply the 3' rule to delins variants, get rid of equals, and substitute
    any slices relative to the observed sequence.
    """
    new_variants = []
    for variant in variants:
        if variant.get("type") == "inversion":
            new_variants.append(copy.deepcopy(variant))
        elif variant.get("type") == "deletion_insertion":
            variant["inserted"] = update_inserted_with_sequences(
                variant["inserted"], sequences
            )
            inserted_sequence = get_inserted_sequence(variant, sequences)
            new_variant = copy.deepcopy(variant)
            shift3 = 0
            shift5 = 0
            if get_location_length(variant["location"]) and not inserted_sequence:
                shift5, shift3 = roll(
                    sequences["reference"],
                    variant["location"]["start"]["position"] + 1,
                    variant["location"]["end"]["position"],
                )
            elif not get_location_length(variant["location"]) and inserted_sequence:
                rolled_sequence = (
                    sequences["reference"][: get_start(variant)]
                    + inserted_sequence
                    + sequences["reference"][get_end(variant) :]
                )
                shift5, shift3 = roll(
                    rolled_sequence,
                    get_start(variant) + 1,
                    get_end(variant) + len(inserted_sequence),
                )
                if shift3:
                    inserted_rolled_sequence = rolled_sequence[
                        get_start(variant)
                        + shift3 : get_end(variant)
                        + shift3
                        + len(inserted_sequence)
                    ]
                    new_variant["inserted"] = [
                        {"sequence": inserted_rolled_sequence, "source": "description"}
                    ]
            shift = shift3 + shift5
            new_variant["location"]["start"]["position"] += shift3
            new_variant["location"]["start"]["shift"] = shift
            new_variant["location"]["end"]["position"] += shift3
            new_variant["location"]["end"]["shift"] = shift
            new_variants.append(new_variant)

    return new_variants


def is_duplication(variant, sequences):
    """
    Note that it works only in the context of the `de_to_hgvs` function flow.
    """

    inserted_sequence = get_inserted_sequence(variant, sequences)
    if len(inserted_sequence) < get_location_length(variant):
        return False
    elif (
        get_start(variant) == get_end(variant)
        and sequences["reference"][
            get_start(variant) - len(inserted_sequence) : get_start(variant)
        ]
        == inserted_sequence
    ):
        return True
    return False


def is_repeat(variant, sequences):
    """
    Note that it works only in the context of the `de_to_hgvs` function flow.
    """
    inserted_sequence = get_inserted_sequence(variant, sequences)
    if len(inserted_sequence) > 2000:
        return False
    repeat_seq, repeat_number = seq_present_before(
        sequences["reference"],
        inserted_sequence,
        get_start(variant["location"]),
        get_end(variant["location"]),
    )
    if repeat_number > 1:
        return True
    return False


def inserted_to_hgvs(inserted):
    new_inserted = []
    for insert in inserted:
        if insert["source"] in ["reference"]:
            new_inserted.append({"source": "reference", "location": insert["location"]})
        else:
            new_inserted.append(
                {"source": "description", "sequence": insert["sequence"]}
            )
    return new_inserted


def delins_to_del(variant):
    return {
        "type": "deletion",
        "source": "reference",
        "location": copy.deepcopy(variant["location"]),
    }


def delins_to_duplication(variant, sequences):
    def _update_shift():
        location = new_variant["location"]
        if location["start"].get("shift"):
            location["start"]["shift"] -= get_location_length(location)
        if location["end"].get("shift"):
            location["end"]["shift"] -= get_location_length(location)

    new_variant = copy.deepcopy(variant)
    inserted_sequence = get_inserted_sequence(variant, sequences)
    new_variant["location"]["start"]["position"] = get_start(
        new_variant["location"]
    ) - len(inserted_sequence)
    new_variant.pop("inserted")
    new_variant["type"] = "duplication"
    _update_shift()
    return new_variant


def delins_to_repeat(variant, sequences):
    new_variant = copy.deepcopy(variant)
    inserted_sequence = get_inserted_sequence(variant, sequences)
    repeat_seq, repeat_number = seq_present_before(
        sequences["reference"],
        inserted_sequence,
        get_start(variant["location"]),
        get_end(variant["location"]),
    )
    shift_left = len(repeat_seq)
    while True:
        if (
            get_start(variant) - len(repeat_seq) > 0
            and sequences["reference"][
                get_start(variant)
                - shift_left
                - len(repeat_seq) : get_start(variant)
                - shift_left
            ]
            == repeat_seq
        ):
            shift_left += len(repeat_seq)
        else:
            break
    repeat_number += shift_left // len(repeat_seq)
    new_variant["location"]["start"]["position"] -= shift_left
    new_variant["type"] = "repeat"
    new_variant["inserted"] = [
        {
            "sequence": repeat_seq,
            "source": "description",
            "repeat_number": {"value": repeat_number, "type": "point"},
        }
    ]
    if new_variant["location"]["start"].get("shift"):
        new_variant["location"]["start"]["shift"] -= shift_left
    if new_variant["location"]["end"].get("shift"):
        new_variant["location"]["end"]["shift"] -= shift_left
    return new_variant


def delins_to_insertion(variant):
    new_variant = copy.deepcopy(variant)
    new_variant["type"] = "insertion"
    return new_variant


def delins_to_substitution(variant, sequences):
    new_variant = copy.deepcopy(variant)
    new_variant["type"] = "substitution"
    new_variant["deleted"] = [
        {
            "sequence": slice_sequence(variant["location"], sequences["reference"]),
            "source": "description",
        }
    ]
    new_variant["inserted"] = inserted_to_hgvs(variant["inserted"])
    return new_variant


def delins_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant["inserted"] = inserted_to_hgvs(variant["inserted"])
    return new_variant


def de_to_hgvs(variants, sequences=None):
    """
    Convert the description extractor variants to an HGVS format (e.g., a
    deletion insertion of one nucleotide is converted to a substitution).
    """
    if len(variants) == 1 and variants[0].get("type") == "equal":
        new_variant = copy.deepcopy(variants[0])
        new_variant.pop("location")
        return [new_variant]

    new_variants = []
    for variant in de_variants_clean(variants, sequences):
        if variant.get("type") == "inversion":
            new_variants.append(copy.deepcopy(variant))
        elif variant.get("type") == "deletion_insertion":
            inserted_sequence = get_inserted_sequence(variant, sequences)
            if len(inserted_sequence) == 0:
                new_variants.append(delins_to_del(variant))
            elif (
                get_location_length(variant["location"]) == len(inserted_sequence) == 1
            ):
                new_variants.append(delins_to_substitution(variant, sequences))
            elif is_repeat(variant, sequences):
                new_variants.append(delins_to_repeat(variant, sequences))
            elif is_duplication(variant, sequences):
                new_variants.append(delins_to_duplication(variant, sequences))
            elif get_start(variant["location"]) == get_end(variant["location"]):
                new_variants.append(delins_to_insertion(variant))
            else:
                new_variants.append(delins_to_delins(variant))

    return new_variants
