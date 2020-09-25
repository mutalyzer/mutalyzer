from .description import Description


def name_check(description_to_normalize):
    description = Description(description_to_normalize)

    description.normalize()

    return description.output()
