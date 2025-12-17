import argparse
from mutalyzer.description_model import get_reference_id, model_to_string
from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer.reference import retrieve_reference, yield_feature_models, is_overlap
from mutalyzer.description_model import get_reference_id, model_to_string
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


def gene_to_transcripts(reference_id: str, gene_symbol: str):
    """Return all transcripts for a given gene in a reference."""
    reference_model, found = get_reference_model(reference_id)
    if not found:
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
    return transcripts

def get_normalized_model(description):
    d = Description(
        description=description,
        only_variants=False,
        sequence=None,
    )
    d.normalize(include_extras=True)

    if d.errors:
        return d
    if not d.references and not d.references.get("observed"):
        d.errors.append({
            "details": "No observed sequence or other error occurred.",
            "source": "input",
        })
    return d

def get_canonical_variants(description):
    """Return the superemals of a given variant description, algebra based extractor."""
    d = get_normalized_model(description)

    if d.errors:
        return d
    observed  = d.references["observed"]["sequence"]["seq"]
    reference = d.references["reference"]["sequence"]["seq"]
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


def convert_description(description, selector_id):
    """Convert a given variant description to a transcript specific description."""
    #TODO:
        # -Discuss if the description should be noramalized or not, now it is normalized first
        # -Add support for transcript descriptison as input (e,g, NM_001127208.3:c.100del)
        #  via chromosomal description?

    # -normalize the description
    # -check if selector is in the affected genes
    d = get_normalized_model(description)
    if d.errors:
        return d

    gene = transcript_to_gene(selector_id)
    if not gene:
        raise ValueError(f"Cannot find annotated gene for transcript {selector_id}.")

    affected_genes = set()
    boundaries = get_canonical_variants_boundaries(description)
    for start, end in boundaries:
        affected_gene = overlap_genes(
            reference_id=get_reference_id(d.corrected_model),
            start=start,
            end=end,
        )
        affected_genes.update(affected_gene)
    if gene not in affected_genes:
        raise ValueError(f"Selector transcript {selector_id} does not belong to affected genes {affected_genes}.")

    converted_descriptions = []
    if d.de_hgvs_internal_indexing_model:
        from_model = d.de_hgvs_internal_indexing_model
    else:
        from_model = None

    try:
        converted_model = to_hgvs_locations(
            model=from_model,
            references=d.references,
            to_coordinate_system="c",
            to_selector_id=selector_id,
            degenerate=True,
        )
    except Exception as e:
        return {"errors": [{"details": str(e)}], "source": "conversion"}
    converted_descriptions.append(model_to_string(converted_model))
    return converted_descriptions


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
    print(convert_description("NC_000004.12(XM_047415843.1):c.100del", "NM_001127208.3"))

    print(overlap_genes(args.reference, args.start, args.end))
    print(gene_to_transcripts(args.reference, "SDHD"))
    # print(overlap_mane_selectors(args.reference, args.start, args.end))