from .description import Description
import time


def normalize(description, only_variants=False, sequence=None):
    t0 = time.time()
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )
    d.normalize()
    t = time.time()
    d.get_chromosomal_descriptions()
    print("get_chromosomal_description:", time.time() - t)
    print("- TOTAL time:", time.time() - t0, "\n\n")
    output = d.output()

    return output
