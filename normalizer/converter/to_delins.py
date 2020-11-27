import copy

from normalizer.util import set_start, get_end, get_location_length


def substitution_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant["type"] = "deletion_insertion"
    return new_variant


def deletion_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant["type"] = "deletion_insertion"
    new_variant["inserted"] = []
    return new_variant


def duplication_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant["type"] = "deletion_insertion"
    new_variant["inserted"] = [
        {"source": "reference", "location": copy.deepcopy(new_variant["location"])}
    ]
    set_start(new_variant["location"], get_end(new_variant["location"]))
    return new_variant


def repeat_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant["type"] = "deletion_insertion"

    # TODO: can this ever be > 1?
    assert len(new_variant["inserted"]) == 1

    ll = get_location_length(new_variant["location"])
    il = len(new_variant["inserted"][0]["sequence"])
    if ll != il:
        raise Exception("Range length and sequence don't match")

    # TODO: it would be nice to match the repeated sequence with the reference,
    #       but I'm not sure if and where that can be done

    new_variant["inserted"][0]["duplication"] = new_variant["inserted"][0]["repeat_number"]["value"]
    return new_variant


def insertion_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant["type"] = "deletion_insertion"
    return new_variant


def inversion_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant["type"] = "deletion_insertion"
    new_variant["inserted"] = [
        {
            "source": "reference",
            "location": copy.deepcopy(new_variant["location"]),
            "inverted": True,
        }
    ]
    return new_variant


def conversion_to_delins(variant):
    new_variant = copy.deepcopy(variant)
    new_variant["type"] = "deletion_insertion"
    return new_variant


def deletion_insertion_to_delins(variant):
    return copy.deepcopy(variant)


def equal_to_delins(variant):
    """
    Only works for variants using internal indexing
    """
    new_variant = copy.deepcopy(variant)
    new_variant["inserted"] = [
        {"source": "reference", "location": copy.deepcopy(new_variant["location"])}
    ]
    new_variant["type"] = "deletion_insertion"
    return new_variant


def variants_to_delins(variants):
    """
    Convert the variants list to its deletion insertion only
    equivalent. It considers that internal indexing is employed.
    """
    new_variants = []
    for variant in variants:
        if variant.get('type') == 'substitution':
            new_variants.append(substitution_to_delins(variant))
        elif variant.get('type') == "deletion":
            new_variants.append(deletion_to_delins(variant))
        elif variant.get('type') == "duplication":
            new_variants.append(duplication_to_delins(variant))
        elif variant.get('type') == "insertion":
            new_variants.append(insertion_to_delins(variant))
        elif variant.get('type') == "inversion":
            new_variants.append(inversion_to_delins(variant))
        elif variant.get('type') == "deletion_insertion":
            new_variants.append(deletion_insertion_to_delins(variant))
        elif variant.get('type') == "equal":
            new_variants.append(equal_to_delins(variant))
        elif variant.get('type') == "repeat":
            new_variants.append(repeat_to_delins(variant))
        else:
            # TODO: Add error.
            print("no supported variant type")

    return new_variants


def to_delins(model):
    new_model = copy.deepcopy(model)
    if new_model.get("variants"):
        new_model["variants"] = variants_to_delins(model["variants"])
    return new_model