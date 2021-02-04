from .description import Description


def normalize(description_to_normalize):
    description = Description(description=description_to_normalize,
                              stop_on_error=True)

    description.normalize()

    return description.normalized_description
