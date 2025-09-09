from .description import Description, view_delins


def view_variants_normalized(d, left=15, right=15, invert=True):
    output = view_delins(
        d.delins_model["variants"],
        d.corrected_model["variants"][::-1] if d.is_inverted() else d.corrected_model["variants"],
        d.get_sequences(),
        left,
        right,
        invert and d.is_inverted(),
    )
    if d.infos:
        output["infos"] = d.infos
    return output


def view_variants(
    description, only_variants=False, sequence=None, left=15, right=15, invert=True
):
    d = Description(
        description=description, only_variants=only_variants, sequence=sequence
    )
    d.to_delins()

    if d.errors:
        return d.output()
    else:
        return view_variants_normalized(d, left, right, invert)
