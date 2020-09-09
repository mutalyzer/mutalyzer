import copy

from normalizer.util import get_end, set_start


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
        if variant['type'] == 'substitution':
            new_variants.append(substitution_to_delins(variant))
        elif variant['type'] == "deletion":
            new_variants.append(deletion_to_delins(variant))
        elif variant['type'] == "duplication":
            new_variants.append(duplication_to_delins(variant))
        elif variant['type'] == "insertion":
            new_variants.append(insertion_to_delins(variant))
        elif variant['type'] == "inversion":
            new_variants.append(inversion_to_delins(variant))
        elif variant['type'] == "deletion_insertion":
            new_variants.append(deletion_insertion_to_delins(variant))
        elif variant['type'] == "equal":
            new_variants.append(equal_to_delins(variant))
        else:
            # TODO: Add error.
            print("no supported variant type")

    return new_variants


def to_delins(model):
    new_model = copy.deepcopy(model)
    if new_model.get("variants"):
        new_model["variants"] = variants_to_delins(model["variants"])
    return new_model