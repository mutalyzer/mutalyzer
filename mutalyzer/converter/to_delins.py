import copy

from mutalyzer.util import get_end, get_location_length, set_start


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
    if variant.get("location") is None:
        return
    new_variant = copy.deepcopy(variant)
    new_variant["inserted"] = [
        {"location": new_variant["location"], "source": "reference"}
    ]
    new_variant["type"] = "deletion_insertion"
    return new_variant


def variant_to_delins(variant):
    """
    Convert the variant to its deletion insertion only equivalent.
    It considers that internal indexing is employed.
    """
    if variant.get("type") == "substitution":
        return substitution_to_delins(variant)
    elif variant.get("type") == "deletion":
        return deletion_to_delins(variant)
    elif variant.get("type") == "duplication":
        return duplication_to_delins(variant)
    elif variant.get("type") == "insertion":
        return insertion_to_delins(variant)
    elif variant.get("type") == "inversion":
        return inversion_to_delins(variant)
    elif variant.get("type") == "conversion":
        return conversion_to_delins(variant)
    elif variant.get("type") == "deletion_insertion":
        return deletion_insertion_to_delins(variant)
    elif variant.get("type") == "equal":
        return equal_to_delins(variant)
    elif variant.get("type") == "repeat":
        return repeat_to_delins(variant)
    # TODO: Add error?
    print("No variant or not supported variant type.")


def variants_to_delins(variants):
    """
    Convert the variants list to its deletion insertion only
    equivalent. It considers that internal indexing is employed.
    """
    new_variants = []
    for variant in variants:
        new_variant = variant_to_delins(variant)
        if new_variant:
            new_variants.append(new_variant)

    return new_variants


def to_delins(model):
    new_model = copy.deepcopy(model)
    if new_model.get("variants"):
        new_model["variants"] = variants_to_delins(model["variants"])
    return new_model
