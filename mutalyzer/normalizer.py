import time

from .description import Description


def normalize(description, only_variants=False, sequence=None):
    # t0 = time.time()
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )
    d.normalize()
    t = time.time()
    d.get_chromosomal_descriptions()
    # print("get_chromosomal_description:", time.time() - t)
    # print("- TOTAL time:", time.time() - t0, "\n\n")
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
        output["delins_model"] = d.delins_model
    return output
