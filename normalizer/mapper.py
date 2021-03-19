from copy import deepcopy

import extractor
from mutalyzer_crossmapper import NonCoding
from mutalyzer_mutator.util import reverse_complement
from mutalyzer_retriever.retriever import NoReferenceError, NoReferenceRetrieved

from .converter import de_to_hgvs
from .converter.to_hgvs_coordinates import (
    crossmap_to_hgvs_setup,
    locations_to_hgvs_locations,
    to_hgvs_locations,
)
from .converter.to_hgvs_indexing import variants_to_internal_indexing
from .description import Description
from .description_model import model_to_string, variants_to_description
from .errors import no_selector_found, reference_not_retrieved
from .reference import (
    get_coordinate_system_from_reference,
    get_coordinate_system_from_selector_id,
    get_reference_model,
    get_reference_mol_type,
    get_selector_model,
)


def _extract_seq(seq, slices):
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


def map_description(description, reference_id, selector_id=None, clean=False):
    d = Description(description)
    d.normalize()
    if d.errors:
        return {"errors": d.errors}
    if d.references and d.references.get("observed"):
        try:
            r_model = get_reference_model(reference_id)
        except (NoReferenceError, NoReferenceRetrieved):
            return {"errors": [reference_not_retrieved(reference_id, [])]}

        if selector_id:
            s_model = get_selector_model(r_model["annotations"], selector_id, True)
            if s_model is None:
                return {"errors": [no_selector_found(reference_id, selector_id, [])]}
            ref_seq = _extract_seq(r_model["sequence"]["seq"], s_model["exon"])
            if s_model["inverted"]:
                ref_seq = reverse_complement(ref_seq)
        else:
            ref_seq = r_model["sequence"]["seq"]
        obs_seq = d.references["observed"]["sequence"]["seq"]

        de_hgvs_variants = de_to_hgvs(
            extractor.describe_dna(ref_seq, obs_seq),
            {"reference": ref_seq, "observed": obs_seq},
        )

        if clean:
            raw_de_hgvs_variants = de_to_hgvs(
                extractor.describe_dna(
                    ref_seq, d.references["reference"]["sequence"]["seq"]
                ),
                {
                    "reference": ref_seq,
                    "observed": d.references["reference"]["sequence"]["seq"],
                },
            )
            de_hgvs_variants = [
                v for v in de_hgvs_variants if v not in raw_de_hgvs_variants
            ]

        if get_reference_mol_type(r_model) in ["mRNA"]:
            return model_to_string(
                to_hgvs_locations(
                    model={
                        "reference": {"id": reference_id},
                        "variants": de_hgvs_variants,
                    },
                    references={"reference": r_model, reference_id: r_model},
                    to_coordinate_system="c",
                    to_selector_id=reference_id,
                    degenerate=True,
                )
            )
        elif selector_id:
            s_model = _convert_selector_locations(s_model)
            return model_to_string(
                to_hgvs_locations(
                    model={
                        "reference": {
                            "id": reference_id,
                            "selector": {"id": selector_id},
                        },
                        "variants": de_hgvs_variants,
                    },
                    to_coordinate_system="c",
                    references={"reference": r_model, reference_id: r_model},
                    degenerate=True,
                    selector_model=s_model,
                )
            )
        elif get_coordinate_system_from_reference(r_model) == "n":
            return model_to_string(
                to_hgvs_locations(
                    model={
                        "reference": {"id": reference_id},
                        "variants": de_hgvs_variants,
                    },
                    references={"reference": r_model, reference_id: r_model},
                    to_selector_id=reference_id,
                    degenerate=True,
                )
            )
        else:
            hgvs_indexing_variants = variants_to_internal_indexing(de_hgvs_variants)
            c_s = get_coordinate_system_from_selector_id(r_model, reference_id)
            crossmap = crossmap_to_hgvs_setup(c_s)
            hgvs_variants = locations_to_hgvs_locations(
                {"variants": hgvs_indexing_variants}, crossmap
            )
            return model_to_string(
                {
                    "reference": {"id": reference_id},
                    "coordinate_system": c_s,
                    "variants": hgvs_variants["variants"],
                }
            )
    else:
        return {"errors": [{"details": "No observed sequence or other error occured."}]}
