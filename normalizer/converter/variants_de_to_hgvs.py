import copy

from normalizer.util import get_end, get_location_length, get_start, roll


def is_deletion(delins_variant):
    if (delins_variant.get("inserted") is None) or (
        len(delins_variant.get("inserted")) == 0
    ):
        return True
    for inserted in delins_variant["inserted"]:
        if get_location_length(inserted["location"]) > 0:
            return False
    return True


def slice_sequence(location, sequence):
    return sequence[get_start(location) : get_end(location)]


def update_inserted_with_sequences(inserted, sequences):
    for insert in inserted:
        if insert["source"] == "observed":
            insert["sequence"] = sequences["observed"][
                get_start(insert["location"]) : get_end(insert["location"])
            ]


def get_inserted_sequence(variant, sequences):
    seq = ""
    if variant.get("inserted"):
        for insert in variant["inserted"]:
            seq += slice_sequence(insert["location"], sequences[insert["source"]])
        return seq
    return seq


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
            update_inserted_with_sequences(variant["inserted"], sequences)
            inserted_sequence = get_inserted_sequence(variant, sequences)
            new_variant = copy.deepcopy(variant)
            shift3 = 0
            if get_location_length(variant["location"]) and not inserted_sequence:
                shift5, shift3 = roll(
                    sequences["reference"],
                    variant["location"]["start"]["position"] + 1,
                    variant["location"]["end"]["position"],
                )
            elif not get_location_length(variant["location"]) and inserted_sequence:
                shift5, shift3 = roll(
                    sequences["reference"][: get_start(variant)]
                    + inserted_sequence
                    + sequences["reference"][get_end(variant) :],
                    get_start(variant) + 1,
                    get_end(variant) + len(inserted_sequence),
                )
            new_variant["location"]["start"]["position"] += shift3
            new_variant["location"]["end"]["position"] += shift3
            new_variants.append(new_variant)

    return new_variants


def is_duplication(variant, sequences):
    """
    Note that it works only in the context of the `de_to_hgvs` function flow.
    """
    inserted_sequence = get_inserted_sequence(variant, sequences)
    if len(inserted_sequence) > get_start(variant):
        return False
    elif (
        sequences["reference"][
            get_start(variant) - len(inserted_sequence) : get_start(variant)
        ]
        == inserted_sequence
    ):
        return True
    return False


def delins_to_del(variant):
    return {
        "type": "deletion",
        "source": "reference",
        "location": copy.deepcopy(variant["location"]),
    }


def delins_to_duplication(variant, sequences):
    new_variant = copy.deepcopy(variant)
    inserted_sequence = get_inserted_sequence(variant, sequences)
    new_variant["location"]["start"]["position"] = get_start(
        new_variant["location"]
    ) - len(inserted_sequence)
    new_variant["type"] = "duplication"
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
            "source": "reference_location",
        }
    ]
    return new_variant


def de_to_hgvs(variants, sequences=None):
    """
    Convert the description extractor variants to an HGVS format (e.g., a
    deletion insertion of one nucleotide is converted to a substitution).
    """
    if len(variants) == 1 and variants[0].get("type") == "equal":
        return [copy.deepcopy(variants[0])]

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
            elif is_duplication(variant, sequences):
                new_variants.append(delins_to_duplication(variant, sequences))
            elif get_start(variant["location"]) == get_end(variant["location"]):
                new_variants.append(delins_to_insertion(variant))
            else:
                new_variants.append(copy.deepcopy(variant))

    return new_variants
