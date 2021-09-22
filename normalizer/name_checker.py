from .description import Description


def name_check(description, only_variants=False, reference_sequence=None):
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=reference_sequence,
    )

    d.normalize()

    return d.output()
