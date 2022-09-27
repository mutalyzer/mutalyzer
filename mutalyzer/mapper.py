import bisect
from copy import deepcopy

import extractor
from mutalyzer_crossmapper import NonCoding
from mutalyzer_mutator import mutate
from mutalyzer_mutator.util import reverse_complement
from mutalyzer_retriever.retriever import NoReferenceError, NoReferenceRetrieved

import mutalyzer.errors as errors

from .converter import de_to_hgvs
from .converter.to_hgvs_coordinates import to_hgvs_locations
from .converter.to_rna import get_location_type
from .description import Description
from .description_model import model_to_string, variants_to_description
from .reference import (
    extract_feature_model,
    get_coordinate_system_from_reference,
    get_coordinate_system_from_selector_id,
    get_internal_selector_model,
    get_reference_model,
    get_selector_feature,
    retrieve_reference,
)
from .util import get_end, get_inserted_sequence, get_start, set_end, set_start


def _slice_seq(seq, slices):
    output = ""
    for s in slices:
        output += seq[s[0] : s[1]]
    return output


def _get_gene_locations(r_model):
    g_l = r_model["annotations"]["features"][0]["location"]
    return g_l["start"]["position"], g_l["end"]["position"]


def _yield_locations(r_model, selector_id):
    for feature in get_selector_feature(r_model["annotations"], selector_id)[
        "features"
    ]:
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
    s_model = get_internal_selector_model(new_r_model["annotations"], selector_id, True)

    if slice_to == "transcript":
        ref_seq = _slice_seq(ref_seq, s_model["exon"])
        new_r_model["sequence"] = {"seq": ref_seq}
        x = NonCoding(s_model["exon"]).coordinate_to_noncoding
        exon_end = [exon[1] for exon in s_model["exon"]]
    if slice_to == "gene":
        g_l = _get_gene_locations(new_r_model)
        ref_seq = _slice_seq(ref_seq, [g_l])
        new_r_model["sequence"] = {"seq": ref_seq}
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
        s_model = get_internal_selector_model(r_model, selector_id, True)
        if s_model["inverted"]:
            ref_seq = reverse_complement(ref_seq)
    return ref_seq


def _get_description(de_hgvs_internal_indexing_variants, r_model, selector_id=None):
    reference = {"id": r_model["annotations"]["id"]}
    if selector_id:
        reference["selector"] = {"id": selector_id}
        c_s = get_coordinate_system_from_selector_id(r_model, selector_id)
    else:
        c_s = get_coordinate_system_from_reference(r_model)
        if c_s in ["c", "n"]:
            selector_id = r_model["annotations"]["id"]

    de_hgvs_model = to_hgvs_locations(
        {
            "reference": reference,
            "coordinate_system": "i",
            "variants": de_hgvs_internal_indexing_variants,
        },
        {"reference": r_model, r_model["annotations"]["id"]: r_model},
        c_s,
        selector_id,
        True,
    )
    return model_to_string(de_hgvs_model)


def _extract_hgvs_internal_model(obs_seq, ref_seq):
    de_variants = extractor.describe_dna(ref_seq, obs_seq)

    return de_to_hgvs(
        de_variants,
        {"reference": ref_seq, "observed": obs_seq},
    )


def _filter(variants, ref_seq1, ref_seq2):
    raw_de_variants = extractor.describe_dna(ref_seq1, ref_seq2)
    seq_variants = de_to_hgvs(
        raw_de_variants,
        {"reference": ref_seq1, "observed": ref_seq2},
    )
    return [v for v in variants if v not in seq_variants]


def _slice_from(variants, exons, sequences):
    sliced = []
    unsliced = []
    x = NonCoding(exons).coordinate_to_noncoding
    for i, v in enumerate(variants):
        slice_v = deepcopy(v)
        if v.get("location"):
            if get_location_type(v["location"], exons, 0, 0) in [
                "same exon",
                "exon exon",
            ]:
                set_start(slice_v["location"], x(get_start(slice_v))[0] - 1)
                set_end(
                    slice_v["location"],
                    x(get_end(slice_v))[0] + x(get_end(slice_v))[1] - 1,
                )
                if slice_v.get("inserted"):
                    slice_v["inserted"] = [
                        {
                            "source": "description",
                            "sequence": get_inserted_sequence(slice_v, sequences),
                        }
                    ]
                sliced.append(slice_v)
            else:
                unsliced.append(i)
    return sliced, unsliced


def map_description(
    description,
    reference_id,
    selector_id=None,
    slice_to=None,
    filter=False,
    len_max=100000,
    diff_max=1000,
):
    # Get the observed sequence
    d = Description(description)
    d.normalize()
    if d.errors:
        return {"errors": d.errors}
    if not d.references and not d.references.get("observed"):
        return {"errors": [{"details": "No observed sequence or other error occured."}]}
    obs_seq = d.references["observed"]["sequence"]["seq"]

    to_r_model = retrieve_reference(reference_id)[0]
    if to_r_model is None:
        return {"errors": [errors.reference_not_retrieved(reference_id, [])]}

    ref_seq2 = d.references["reference"]["sequence"]["seq"]

    if d.only_equals() or d.no_operation():
        variants = []
    else:
        variants = d.delins_model["variants"]

    if slice_to == "transcript":
        selector_model = d.get_selector_model()
        if selector_model:
            sliced_variants, unsliced_variants = _slice_from(
                variants,
                selector_model["exon"],
                d.get_sequences(),
            )
            if unsliced_variants:
                errs = []
                for v in unsliced_variants:
                    errs.append(
                        errors.location_slice(
                            d.corrected_model["variants"][v]["location"]
                        )
                    )
                return {"errors": errs}
            from_r_model = _get_reference_model(
                d.references["reference"], d.get_selector_id(), slice_to
            )
            ref_seq2 = from_r_model["sequence"]["seq"]
            obs_seq = mutate({"reference": ref_seq2}, sliced_variants)
    elif slice_to == "gene":
        gene = extract_feature_model(
            d.references["reference"]["annotations"], d.get_selector_id()
        )[0]
        if gene:
            new_r_model = {"annotations": deepcopy(gene)}
            g_l = _get_gene_locations(new_r_model)
            ref_seq2 = _slice_seq(d.references["reference"]["sequence"]["seq"], [g_l])
            sliced_variants, unsliced_variants = _slice_from(
                variants, [g_l], {"reference": ref_seq2}
            )
            if unsliced_variants:
                errs = []
                for v in unsliced_variants:
                    errs.append(
                        errors.location_slice(
                            d.corrected_model["variants"][v]["location"]
                        )
                    )
                return {"errors": errs}
            obs_seq = mutate({"reference": ref_seq2}, sliced_variants)
    elif slice_to is not None:
        return {"errors": [errors.slice_option(slice_to)]}

    if selector_id:
        s_model = get_internal_selector_model(
            to_r_model["annotations"], selector_id, True
        )
        if s_model is None:
            return {"errors": [errors.no_selector_found(reference_id, selector_id, [])]}
        if s_model["inverted"] and not (
            d.get_selector_model() and d.get_selector_model()["inverted"]
        ):
            obs_seq = reverse_complement(obs_seq)
            ref_seq2 = reverse_complement(ref_seq2)
    if slice_to:
        to_r_model = _get_reference_model(to_r_model, selector_id, slice_to)

    ref_seq1 = to_r_model["sequence"]["seq"]

    if len(ref_seq1) > len_max:
        return {"errors": [errors.sequence_length(ref_seq1, len_max)]}
    if len(obs_seq) > len_max:
        return {"errors": [errors.sequence_length(obs_seq, len_max)]}

    if len(ref_seq1) < len(obs_seq) and abs(len(ref_seq1) - len(obs_seq)) > diff_max:
        return {
            "errors": [
                errors.lengths_difference(abs(len(ref_seq1) - len(obs_seq)), diff_max)
            ]
        }
    # Get the description extractor hgvs internal indexing variants
    variants = _extract_hgvs_internal_model(obs_seq, ref_seq1)

    if filter:
        raw_de_variants = extractor.describe_dna(ref_seq1, ref_seq2)
        seq_variants = de_to_hgvs(
            raw_de_variants,
            {"reference": ref_seq1, "observed": ref_seq2},
        )
        if not (len(seq_variants) == 1 and seq_variants[0]["type"] == "equal") and [
            v for v in seq_variants if v not in variants
        ]:
            return {
                "errors": [{"code": "EMAPFILTER", "details": "Unsuccessful filtering."}]
            }
        variants = [v for v in variants if v not in seq_variants]

    return {"mapped_description": _get_description(variants, to_r_model, selector_id)}
