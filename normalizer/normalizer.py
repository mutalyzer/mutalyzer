from .description import Description


def normalize(description_to_normalize):
    description = Description(description_to_normalize)

    description.normalize()

    return description.normalized_description
