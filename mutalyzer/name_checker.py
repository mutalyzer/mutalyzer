from .description import Description


def normalizer_(description, only_variants=False, sequence=None):
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )

    d.normalize()

    return d.output()
