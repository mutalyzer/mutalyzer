from .description import Description


def mutate(description):
    d = Description(description)

    d.normalize()

    if d.references.get("observed"):
        return d.references["observed"]
    else:
        status = d.output()
        output = {}
        if status.get("errors"):
            output["errors"] = status["errors"]
        if status.get("infos"):
            output["infos"] = status["infos"]
        return output


def mutate_sequence(seq, variants):
    pass
