from copy import deepcopy

import extractor
from mutalyzer_crossmapper import NonCoding
from mutalyzer_mutator.util import reverse_complement
from mutalyzer_retriever.retriever import NoReferenceError, NoReferenceRetrieved

from .converter import de_to_hgvs
from .converter.to_hgvs_coordinates import to_hgvs_locations
from .description import Description
from .description_model import model_to_string, variants_to_description
from .errors import no_selector_found, reference_not_retrieved
from .reference import (
    extract_feature_model,
    get_coordinate_system_from_reference,
    get_coordinate_system_from_selector_id,
    get_feature,
    get_reference_model,
    get_selector_model,
)


def _slice_seq(seq, slices):
    output = ""
    for slice in slices:
        output += seq[slice[0] : slice[1]]
    return output


def _convert_selector_locations(s_model):
    output = deepcopy(s_model)
    exon = []
    x = NonCoding(s_model["exon"], s_model["inverted"]).coordinate_to_noncoding
    for e in s_model["exon"]:
        exon.append((x(e[0])[0] - 1, x(e[1])[0]))
    output["exon"] = exon
    output["cds"] = [(x(s_model["cds"][0][0])[0] - 1, x(s_model["cds"][0][1])[0])]
    return output


def _get_gene_locations(r_model):
    g_l = r_model["annotations"]["features"][0]["location"]
    return g_l["start"]["position"], g_l["end"]["position"]


def _yield_locations(r_model, selector_id):
    for feature in get_feature(r_model["annotations"], selector_id)["features"]:
        if feature.get("location"):
            yield feature["location"], feature["type"]


def _get_reference_model(r_model, selector_id=None, slice_to=None):
    """
    Modify the reference model based on the `selector_id`
    and the `slice_to` parameters.

    :param r_model:
    :param selector_id:
    :param slice_to:
    :return:
    """
    if selector_id is None:
        return r_model

    new_r_model = {
        "annotations": deepcopy(
            extract_feature_model(r_model["annotations"], selector_id)[0]
        )
    }
    ref_seq = r_model["sequence"]["seq"]

    s_model = get_selector_model(new_r_model["annotations"], selector_id, True)

    if slice_to == "transcript":
        new_r_model["sequence"] = {
            "seq": _slice_seq(r_model["sequence"]["seq"], s_model["exon"])
        }
        x = NonCoding(s_model["exon"], s_model["inverted"]).coordinate_to_noncoding
        exon_end = [exon[1] for exon in s_model["exon"]]
    if slice_to == "gene":
        g_l = _get_gene_locations(new_r_model)
        new_r_model["sequence"] = {"seq": _slice_seq(r_model["sequence"]["seq"], [g_l])}
        x = NonCoding([g_l]).coordinate_to_noncoding
        exon_end = [g_l[1]]

    for loc, feature_type in _yield_locations(new_r_model, selector_id):
        loc["start"]["position"] = x(loc["start"]["position"])[0] - 1
        if feature_type == "exon" and loc["end"]["position"] in exon_end:
            loc["end"]["position"] = x(loc["end"]["position"])[0]
        else:
            loc["end"]["position"] = x(loc["end"]["position"])[0] - 1
    return new_r_model


def _get_ref_seq(r_model, selector_id=None):
    """

    :param r_model:
    :param selector_id:
    :param slice_to:
    :return:
    """
    ref_seq = r_model["sequence"]["seq"]
    if selector_id:
        s_model = get_selector_model(r_model, selector_id, True)
        if s_model["inverted"]:
            ref_seq = reverse_complement(ref_seq)
    return ref_seq


def _clean(de_variants, ref_seq1, ref_seq2):
    raw_de_variants = extractor.describe_dna(ref_seq1, ref_seq2)
    return [v for v in de_variants if v not in raw_de_variants]


def _extract_description(de_variants, obs_seq, r_model, selector_id=None):
    ref_seq = r_model["sequence"]["seq"]
    reference = {"id": r_model["annotations"]["id"]}
    if selector_id:
        reference["selector"] = {"id": selector_id}
        c_s = get_coordinate_system_from_selector_id(r_model, selector_id)
    else:
        c_s = get_coordinate_system_from_reference(r_model)
        if c_s in ["c", "n"]:
            selector_id = r_model["annotations"]["id"]

    de_hgvs_internal_indexing_model = {
        "reference": reference,
        "coordinate_system": "i",
        "variants": de_to_hgvs(
            de_variants,
            {"reference": ref_seq, "observed": obs_seq},
        ),
    }
    de_hgvs_model = to_hgvs_locations(
        de_hgvs_internal_indexing_model,
        {"reference": r_model, r_model["annotations"]["id"]: r_model},
        c_s,
        selector_id,
        True,
    )
    return model_to_string(de_hgvs_model)


def map_description(
    description,
    reference_id,
    selector_id=None,
    slice_to=None,
    clean=False,
):
    # Get the observed sequence
    d = Description(description)
    d.normalize()
    if d.errors:
        return {"errors": d.errors}
    if not d.references and not d.references.get("observed"):
        return {"errors": [{"details": "No observed sequence or other error occured."}]}
    obs_seq = d.references["observed"]["sequence"]["seq"]

    # Get the reference_model
    try:
        r_model = get_reference_model(reference_id)
    except (NoReferenceError, NoReferenceRetrieved):
        return {"errors": [reference_not_retrieved(reference_id, [])]}

    if selector_id:
        s_model = get_selector_model(r_model["annotations"], selector_id, True)
        if s_model is None:
            return {"errors": [no_selector_found(reference_id, selector_id, [])]}

    r_model = _get_reference_model(r_model, selector_id, slice_to)

    de_variants = extractor.describe_dna(r_model["sequence"]["seq"], obs_seq)
    if clean:
        de_variants = _clean(
            de_variants,
            r_model["sequence"]["seq"],
            d.references["reference"]["sequence"]["seq"],
        )

    return _extract_description(de_variants, obs_seq, r_model, selector_id)
