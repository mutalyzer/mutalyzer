"""Convert from SPDI to HGVS."""

from extractor import describe_dna
from mutalyzer_mutator import mutate
from mutalyzer_retriever.reference import (
    get_assembly_chromosome_accession,
    get_assembly_id,
)
from mutalyzer_retriever.retriever import get_chromosome_from_selector
from mutalyzer_spdi_parser.convert import to_hgvs_internal_model as spdi_to_hgvs

from .converter.to_hgvs_coordinates import to_hgvs_locations
from .converter.variants_de_to_hgvs import de_to_hgvs
from .description_model import model_to_string
from .errors import out_of_boundary_greater, reference_not_retrieved, sequence_mismatch
from .reference import get_coordinate_system_from_reference, retrieve_reference
from .util import get_end, get_start


def _errors_check(variant, ref_seq):
    """
    Check if the deleted sequence matches the reference sequence and
    if the locations are within the boundaries of te reference sequence.

    :param variant: Description model.
    :param ref_seq: Reference sequence
    :return: A list with errors.
    """
    errors = []
    if variant.get("deleted") and variant["deleted"][0].get("sequence"):
        del_ref_seq = ref_seq[get_start(variant) : get_end(variant)]
        del_seq = variant["deleted"][0]["sequence"]
        if del_seq != del_ref_seq:
            errors.append(sequence_mismatch(del_ref_seq, del_seq, ["variants", 0]))
    if get_start(variant) > len(ref_seq):
        errors.append(
            out_of_boundary_greater(
                variant["location"]["start"],
                get_start(variant) - len(ref_seq),
                len(ref_seq),
                ["variants", 0, "location"],
            )
        )
    if get_end(variant) > len(ref_seq):
        errors.append(
            out_of_boundary_greater(
                variant["location"]["start"],
                get_end(variant) - len(ref_seq),
                len(ref_seq),
                ["variants", 0, "location"],
            )
        )
    return errors


def _get_ids(model):
    """
    Correct `assembly(chromosome_number)`, e.g., `GRCh38(chr11)`.
    """
    r_id = model["reference"]["id"]
    s_id = None
    if model["reference"].get("selector") and model["reference"]["selector"].get("id"):
        s_id = model["reference"]["selector"]["id"]

    a_id = get_assembly_id(r_id)
    if a_id:
        chromosome_id = get_assembly_chromosome_accession(r_id, s_id)
        if chromosome_id:
            r_id = chromosome_id
            s_id = None
        else:
            chromosome_id = get_chromosome_from_selector(a_id, s_id)
            if chromosome_id:
                r_id = chromosome_id
                s_id = None
    return r_id, s_id


def spdi_converter(description):
    """
    Convert an SPDI description to a normalized HGVS description.

    :param description: An SPDI description.
    :return: The normalized HGVS description and the input model.
    """
    try:
        model = spdi_to_hgvs(description)
    except Exception as e:
        if "Unexpected end-of-input" in str(e):
            return {"errors": [{"code": "ESYNTAXUEOF", "details": "Unexpected end of input."}]}
        if "No terminal matches" in str(e):
            return {"errors": [{"code": "ESYNTAXUC", "details": "Unexpected character."}]}
        return {"errors": [{"code": "EUNKNOWN"}]}

    r_id, s_id = _get_ids(model)
    r_m = retrieve_reference(r_id, s_id)[0]

    if r_m is None:
        return {"errors": [reference_not_retrieved(r_id, ["reference", "id"])]}

    errors = _errors_check(model["variants"][0], r_m["sequence"]["seq"])
    if errors:
        return {"errors": errors}

    obs_seq = mutate({"reference": r_m["sequence"]["seq"]}, model["variants"])

    d_v = describe_dna(r_m["sequence"]["seq"], obs_seq)
    d_h_m = {
        "variants": de_to_hgvs(
            d_v, {"reference": r_m["sequence"]["seq"], "observed": obs_seq}
        ),
        "reference": {"id": r_id},
        "coordinate_system": "i",
    }

    c_s = get_coordinate_system_from_reference(r_m)
    if c_s in ["c", "n"]:
        n_m = to_hgvs_locations(d_h_m, {r_id: r_m}, c_s, r_id)
    else:
        n_m = to_hgvs_locations(d_h_m, {r_id: r_m}, c_s)
    n_d = model_to_string(n_m)

    return {"input_model": model, "normalized_description": n_d}
