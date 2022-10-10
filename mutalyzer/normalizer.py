from .description import Description


def normalize(description, only_variants=False, sequence=None):
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )

    d.normalize()

    return d.output()
