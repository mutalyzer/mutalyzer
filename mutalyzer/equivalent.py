import argparse
from mutalyzer_hgvs_parser import to_model
from mutalyzer.description_model import get_reference_id, get_selector_id, model_to_string
from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer.reference import retrieve_reference, yield_feature_models, is_overlap
from mutalyzer.description_model import get_reference_id, model_to_string
from mutalyzer.position_converter import position_convert
from mutalyzer.description import Description
from algebra.lcs.lcs_graph import LCSgraph
from algebra.extractor.extractor import canonical


import pprint


def is_mane(feature_model):
    """Check if a feature model is a MANE Select transcript model."""
    if (
        feature_model.get("qualifiers")
        and feature_model["qualifiers"].get("tag")
        and "MANE" in feature_model["qualifiers"]["tag"]
    )   :
        return True


def get_reference_model(reference_id):
    """Retrieve reference model using reference id."""
    reference_model, _ = retrieve_reference(reference_id)
    return reference_model, reference_model is not None


def overlap_mane_selectors(reference_id, start, end):
    """Return overlapping MANE Select transcripts for a given reference and position."""
    reference_model, found = get_reference_model(reference_id)
    if not found:
        return list()
    mane_selectors = list()
    reference_annotations = reference_model.get("annotations", {})
    for feature in yield_feature_models(reference_annotations):
        if "rna" in feature.get("type").lower() and is_overlap(feature, start, end) and is_mane(feature):
            mane_selectors.append(feature.get("id"))
    return mane_selectors


def annotated_genes(reference_id):
    """Return all annotated genes for a given reference."""
    reference_model, found = get_reference_model(reference_id)
    if not found:
        return list()
    annotated_genes = list()
    reference_annotations = reference_model.get("annotations", {})
    for feature in yield_feature_models(reference_annotations):
        if "gene" in feature.get("type").lower():
            annotated_genes.append(feature.get("id"))
    return annotated_genes

def overlap_genes(reference_id, start, end):
    """"Return overlapping genes for a given reference and position."""
    # retrieve reference model using reference id
    # loop over the features of the model and check if there are features of interest overlapping the position
    #   - feature type: check if is gene
    #   - location range: contains input
    # retrun: a list of overlap genes

    reference_model, found = get_reference_model(reference_id)
    if not found:
        return list()
    overlap_genes = list()
    reference_annotations = reference_model.get("annotations", {})
    for feature in yield_feature_models(reference_annotations):
        if "gene" in feature.get("type").lower() and is_overlap(feature, start, end):
            overlap_genes.append(feature.get("id"))
    return overlap_genes


def annotated_transcripts(reference_id: str, gene_symbol: str):
    """Return all annotated transcripts for a given gene in a reference."""
    reference_model, found = get_reference_model(reference_id)
    if not found or not gene_symbol:
        return list()
    transcripts = list()
    reference_annotations = reference_model.get("annotations", {})
    for feature in yield_feature_models(reference_annotations):
        if (
            "gene" == feature.get("type")
            and gene_symbol.lower() == feature.get("id",[]).lower()
        ):
            for sub_feature in yield_feature_models(feature):
                if "rna" in sub_feature.get("type").lower():
                    transcripts.append(sub_feature.get("id"))
            break
    return transcripts


def get_normalized_model(description):
    d = Description(
        description=description,
        only_variants=False,
        sequence=None,
    )
    d.normalize(include_extras=True)
    return d

def get_canonical_variants(description):
    """Return the superemals of a given variant description, algebra based extractor."""
    normalized_m = get_normalized_model(description)

    if normalized_m.errors:
        return {"errors": normalized_m.errors, "source": "input"}
    if not normalized_m.references and not normalized_m.references.get("observed"):
        return {
            "errors": [{"details": "No observed sequence or other error occurred."}],
            "source": "input",
        }
    observed  = normalized_m.references["observed"]["sequence"]["seq"]
    reference = normalized_m.references["reference"]["sequence"]["seq"]
    Graph = LCSgraph.from_sequence(reference, observed)
    return canonical(Graph)


def get_canonical_variants_string(description):
    variants = get_canonical_variants(description)
    return [repr(v) for v in variants]


def get_canonical_variants_boundaries(description):
    variants = get_canonical_variants(description)
    boundaries = []
    for variant in variants:
        boundaries.append((variant.start, variant.end))
    return boundaries


def transcript_to_gene(transcript_id: str):
    """Return the gene symbol for a given transcript id."""
    # retrieve reference model
    # extract gene feature, return gene symbol
    transcript_model, found = get_reference_model(transcript_id)
    if not found:
        return None
    transcript_annotations = transcript_model.get("annotations", {})
    for feature in yield_feature_models(transcript_annotations):
        if "gene" in feature.get("type").lower():
            return feature.get("id")
    return None


def convert_to_selector_description(description, selector_id):
    """Convert a given variant description to a transcript specific description."""
    #DISCUSSIONS:
        # - Should we normalize the input or the output description?
        #   - Current implementation normalizes the output description

    d_model = to_model(description=description)
    p_c = position_convert(
        description_model=d_model, to_selector_id=selector_id, include_overlapping=False
    )
    if p_c.get("errors") or p_c.get("infos"):
        return p_c
    converted_d = model_to_string(p_c["converted_model"])
    normalized_m = get_normalized_model(converted_d)

    if normalized_m.normalized_description:
        return normalized_m.normalized_description


def convert_to_genomic_description(description):
    """Convert a given variant description to a genomic specific description."""
    d_model = to_model(description=description)
    reference_id = get_reference_id(d_model)
    selector_id = get_selector_id(d_model)
    if not selector_id:
        return {
            "errors": [{"details": "No selector id found in the input description."}],
            "source": "input",
        }
    if selector_id == reference_id:
        return {
            "errors": [{"details": "The selector id is the same as the reference id in the input description."}],
            "source": "input",
        }

    p_c = position_convert(
        description_model=d_model, to_coordinate_system="g", include_overlapping=False
    )
    if p_c.get("errors") or p_c.get("infos"):
        return p_c
    converted_d = model_to_string(p_c["converted_model"])
    normalized_m = get_normalized_model(converted_d)

    if normalized_m.errors:
        return {"errors": normalized_m.errors, "source": "input"}
    elif normalized_m.normalized_description:
        return normalized_m.normalized_description


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate equivalent variant descriptions for a transcript."
    )

    parser.add_argument(
        "reference",
        help="A reference sequence accession (e.g., 'NC_000023.11')"
    )
    parser.add_argument(
        "start",
        type=int,
        help="A location point at the start of the mutation"
    )
    parser.add_argument(
        "end",
        type=int,
        help="A location point at the end of the mutation"
    )


    args = parser.parse_args()
    print(convert_to_selector_description("NC_000004.12(XM_047415843.1):c.100del", "NM_001127208.3"))
    print(convert_to_genomic_description("NM_001127208.3(NM_001127208.3):c.100del"))
    print(convert_to_genomic_description("NC_000004.12(XM_047415843.1):c.100del"))
    print(annotated_genes(args.reference))
    print(overlap_genes(args.reference, args.start, args.end))
    print(annotated_transcripts("NC_000011.10", "SDHD"))
    print(get_canonical_variants("NC_000011.10(NM_003002.4):c.100del"))
    # print(overlap_mane_selectors(args.reference, args.start, args.end))