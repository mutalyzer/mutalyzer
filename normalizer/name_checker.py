from .description import Description


def name_check(description):
    d = Description(description)

    d.normalize()

    return d.output()
