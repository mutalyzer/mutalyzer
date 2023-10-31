"""Interfaces to obtain the normalized (canonical) variant representations or
the delins model of an input description."""

from .description import Description
from .util import construct_sequence


def normalize(description, only_variants=False, sequence=None):
    """
    Obtain the normalized (canonical) variant representation.
    """
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )
    d.normalize()
    d.get_chromosomal_descriptions()
    output = d.output()
    return output


def delins_model(description, only_variants=False, sequence=None):
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )
    d.to_delins()
    output = d.output()
    if d.delins_model:
        for variant in d.delins_model["variants"]:
            if variant.get("inserted"):
                for inserted in variant.get("inserted"):
                    if not inserted.get("sequence"):
                        inserted["sequence"] = construct_sequence(
                            [inserted], d.get_sequences()
                        )
        output["delins_model"] = d.delins_model
    return output
