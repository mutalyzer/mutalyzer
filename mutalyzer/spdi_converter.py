from extractor import describe_dna
from mutalyzer_mutator import mutate
from mutalyzer_spdi_parser.convert import to_hgvs_internal_model as spdi_to_hgvs

from .converter.to_hgvs_coordinates import to_hgvs_locations
from .converter.variants_de_to_hgvs import de_to_hgvs
from .description_model import model_to_string
from .reference import get_coordinate_system_from_reference, retrieve_reference


def spdi_converter(description=None, description_model=None):
    if description:
        m = spdi_to_hgvs(description)
    elif description_model:
        m = description_model

    # TODO: Move in the spdi-parser or support for non inserted in mutator.
    if not m["variants"][0].get("inserted"):
        m["variants"][0]["inserted"] = []

    r_m = retrieve_reference(m["reference"]["id"])[0]
    c_s = get_coordinate_system_from_reference(r_m)

    ref_seq = r_m["sequence"]["seq"]
    obs_seq = mutate({"reference": ref_seq}, m["variants"])
    d_v = describe_dna(r_m["sequence"]["seq"], obs_seq)
    d_h_m = {
        "variants": de_to_hgvs(
            d_v, {"reference": r_m["sequence"]["seq"], "observed": obs_seq}
        ),
        "reference": {"id": m["reference"]["id"]},
        "coordinate_system": "i",
    }
    if c_s == "c":
        n_m = to_hgvs_locations(
            d_h_m, {m["reference"]["id"]: r_m}, c_s, m["reference"]["id"]
        )
    else:
        n_m = to_hgvs_locations(d_h_m, {m["reference"]["id"]: r_m}, c_s)
    n_d = model_to_string(n_m)

    return {"input_model": m, "normalized_description": n_d}
