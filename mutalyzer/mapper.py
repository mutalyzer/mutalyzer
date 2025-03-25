from copy import deepcopy

import extractor
from mutalyzer_mutator import mutate
from mutalyzer_mutator.util import reverse_complement
from mutalyzer_retriever.reference import (
    get_assembly_chromosome_accession,
    get_assembly_id,
)
from mutalyzer_retriever.retriever import extract_feature_model

from mutalyzer import errors

from .converter import de_to_hgvs
from .converter.extras import (
    convert_reference_model,
    convert_to_exons,
    get_gene_locations,
    get_mane_tag,
)
from .converter.to_hgvs_coordinates import to_hgvs_locations
from .description import Description
from .description_model import model_to_string
from .reference import (
    get_coordinate_system_from_reference,
    get_coordinate_system_from_selector_id,
    get_internal_selector_model,
    get_only_selector_id,
    retrieve_reference,
)
from .util import slice_seq


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


def resolve_reference_id(reference_id, selector_id, slice_to, assembly_id, d):
    """
    Resolves `reference_id` based on the annotations retrieved from the input description
    and updates `selector_id` and `slice_to` if necessary.
    """
    if reference_id is None or assembly_id:
        model = d.references["reference"]
        qualifiers = model.get("annotations", {}).get("qualifiers", {})
        chromosome = qualifiers.get("chromosome")
        if chromosome:
            assembly = assembly_id or "GRCH38"
            chr_id = get_assembly_chromosome_accession(assembly, f"chr{chromosome}")
            if chr_id:
                reference_id = chr_id
                if selector_id is None:
                    selector_id = d.get_selector_id()
                if slice_to is None:
                    slice_to = "transcript"
    return reference_id, selector_id, slice_to


def map_description(
    description,
    reference_id=None,
    selector_id=None,
    slice_to=None,
    filter_out=False,
    len_max=110000,
    diff_max=1000,
):
    d = Description(description)
    d.normalize(include_extras=False)
    if d.errors:
        return {"errors": d.errors, "source": "input"}
    if not d.references and not d.references.get("observed"):
        return {
            "errors": [{"details": "No observed sequence or other error occurred."}],
            "source": "input",
        }

    obs_seq = d.references["observed"]["sequence"]["seq"]

    assembly_id = get_assembly_id(reference_id)
    reference_id, selector_id, slice_to = resolve_reference_id(
        reference_id, selector_id, slice_to, assembly_id, d
    )
    to_r_model = retrieve_reference(reference_id, selector_id)[0]

    if to_r_model is None:
        return {
            "errors": [errors.reference_not_retrieved(reference_id, [])],
            "source": "input",
        }

    ref_seq_from = d.references["reference"]["sequence"]["seq"]

    if d.only_equals() or d.no_operation():
        variants = []
    else:
        variants = d.delins_model["variants"]

    if slice_to == "transcript":
        selector_model = d.get_selector_model()
        if selector_model:
            exons = selector_model["exon"]
            if all(exons[i][1] == exons[i + 1][0] for i in range(len(exons) - 1)):
                converted_variants, skipped_variants = variants, []
            else:
                converted_variants, skipped_variants = convert_to_exons(
                    variants,
                    selector_model["exon"],
                    d.get_sequences(),
                )
            if skipped_variants:
                return {
                    "errors": [errors.location_slice(d.corrected_model["variants"][v]["location"])for v in skipped_variants],
                    "source": "input",
                }

            from_r_model = convert_reference_model(
                d.references["reference"], d.get_selector_id(), slice_to
            )
            ref_seq_from = from_r_model["sequence"]["seq"]
            obs_seq = mutate({"reference": ref_seq_from}, converted_variants)
        if (
            selector_id is None
            and get_coordinate_system_from_reference(to_r_model) == "c"
            and get_only_selector_id(to_r_model) == reference_id
        ):
            selector_id = reference_id
    elif slice_to == "gene":
        gene = extract_feature_model(
            d.references["reference"]["annotations"], d.get_selector_id()
        )[0]
        if gene:
            new_r_model = {"annotations": deepcopy(gene)}
            g_l = get_gene_locations(new_r_model)
            ref_seq_from = slice_seq(
                d.references["reference"]["sequence"]["seq"], [g_l]
            )
            converted_variants, skipped_variants = convert_to_exons(
                variants, [g_l], {"reference": ref_seq_from}
            )
            if skipped_variants:
                return {
                    "errors": [errors.location_slice(d.corrected_model["variants"][v]["location"]) for v in skipped_variants],
                    "source": "input",
                }

            obs_seq = mutate({"reference": ref_seq_from}, converted_variants)
    elif slice_to is not None:
        return {"errors": [errors.slice_option(slice_to)], "source": "input"}

    if selector_id:
        s_model = get_internal_selector_model(
            to_r_model["annotations"], selector_id, True
        )
        if s_model is None:
            return {
                "errors": [errors.no_selector_found(reference_id, selector_id, [])],
                "source": "input",
            }
        if d.get_selector_model() and (
            s_model["inverted"] ^ d.get_selector_model()["inverted"]
        ):
            obs_seq = reverse_complement(obs_seq)
            ref_seq_from = reverse_complement(ref_seq_from)
    if slice_to:
        to_r_model = convert_reference_model(to_r_model, selector_id, slice_to)

    ref_seq_to = to_r_model["sequence"]["seq"]

    ref_len = len(ref_seq_to)
    obs_len = len(obs_seq)
    if len_max is not None:
        if ref_len > len_max:
            return {
                "errors": [errors.sequence_length(ref_seq_to, len_max)],
                "source": "input",
            }
        if obs_len > len_max:
            return {
                "errors": [errors.sequence_length(obs_seq, len_max)],
                "source": "input",
            }
    if diff_max is not None:
        length_diff = abs(ref_len - obs_len)
        if ref_len < obs_len and length_diff > diff_max:
            return {
                "errors": [errors.lengths_difference(length_diff, diff_max)],
                "source": "input",
            }

    # Get the description extractor hgvs internal indexing variants
    variants = _extract_hgvs_internal_model(obs_seq, ref_seq_to)

    filtered_variants = False
    unfiltered_mapped_description = None
    reference_sequences_description = None
    if filter_out:
        raw_de_variants = extractor.describe_dna(ref_seq_to, ref_seq_from)
        seq_variants = de_to_hgvs(
            raw_de_variants,
            {"reference": ref_seq_to, "observed": ref_seq_from},
        )
        unfiltered_mapped_description = _get_description(variants, to_r_model, selector_id)
        reference_sequences_description = _get_description(seq_variants, to_r_model, selector_id)
        variants = [v for v in variants if v not in seq_variants]
        if variants != seq_variants:
            filtered_variants = True

    mapped_description = _get_description(variants, to_r_model, selector_id)
    m_d = Description(mapped_description)
    m_d.normalize(include_extras=False)
    m_d.construct_genomic_equivalent()
    if m_d.errors:
        return {"errors": m_d.errors, "source": "output"}
    output = {"mapped_description": m_d.normalized_description}
    if m_d.get_selector_model() and m_d.is_selector_model_valid():
        tag = get_mane_tag(m_d.get_selector_model())
        if tag:
            output["tag"] = tag
    if (
        m_d.equivalent
        and m_d.equivalent.get("g")
        and len(m_d.equivalent["g"]) == 1
        and m_d.equivalent["g"][0].get("description")
    ):
        output["genomic_description"] = m_d.equivalent["g"][0]["description"]
    output["ref_seq_differences"] = ref_seq_to != ref_seq_from
    if filter_out:
        output["filtered_variants"] = filtered_variants
    if unfiltered_mapped_description:
        output["unfiltered_mapped_description"] = unfiltered_mapped_description
    if reference_sequences_description:
        output["reference_sequences_description"] = reference_sequences_description

    return output
