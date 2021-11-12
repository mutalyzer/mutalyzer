from mutalyzer_hgvs_parser import to_model
from mutalyzer_mutator import mutate as mutalyzer_mutator

from .converter.to_delins import to_delins
from .converter.to_internal_coordinates import to_internal_coordinates
from .converter.to_internal_indexing import to_internal_indexing
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
    variants_model = {
        "coordinate_system": "g",
        "variants": to_model(variants, "variants"),
    }
    sequence = {"sequence": {"seq": seq}}
    delins_model = to_delins(
        to_internal_indexing(to_internal_coordinates(variants_model, sequence))
    )
    return mutalyzer_mutator({"reference": seq}, delins_model["variants"])
