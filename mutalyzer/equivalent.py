import argparse
from mutalyzer.description_model import get_reference_id, model_to_string
from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer.reference import retrieve_reference, yield_feature_models, is_overlap
from mutalyzer.description_model import get_reference_id, model_to_string
from mutalyzer.description import Description
from algebra.lcs.lcs_graph import LCSgraph
from algebra.extractor.extractor import canonical


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
        return set()
    mane_selectors = set()
    reference_annotations = reference_model.get("annotations", {})
    for feature in yield_feature_models(reference_annotations):
        if "rna" in feature.get("type").lower() and is_overlap(feature, start, end) and is_mane(feature):
            mane_selectors.add(feature.get("id"))
    return mane_selectors


def annotated_transcripts(reference_id):
    """Return all annotated transcripts for a given reference."""
    reference_model, found = get_reference_model(reference_id)
    if not found:
        return set()
    annotated_transcripts = set()
    reference_annotations = reference_model.get("annotations", {})
    for feature in yield_feature_models(reference_annotations):
        if "rna" in feature.get("type").lower():
            annotated_transcripts.add(feature.get("id"))
    return annotated_transcripts

def overlap_transcripts(reference_id, start, end):
    """"Return overlapping transcripts for a given reference and position."""
    # retrieve reference model using reference id
    # loop over the features of the model and check if there are features of interest overlapping the position
    #   - feature type: check if is transcript
    #   - location range: contains input
    # retrun: a set of overlap transcripts

    reference_model, found = get_reference_model(reference_id)
    if not found:
        return set()
    overlap_transcripts = set()
    reference_annotations = reference_model.get("annotations", {})
    for feature in yield_feature_models(reference_annotations):
        if "rna" in feature.get("type").lower() and is_overlap(feature, start, end):
            overlap_transcripts.add(feature.get("id"))
    return overlap_transcripts


def get_normalized_model(description):
    d = Description(
        description=description,
        only_variants=False,
        sequence=None,
    )
    d.normalize(include_extras=False)

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


def convert_description(description, selectors):
    """Convert a given variant description to a transcript specific description."""
    d = get_normalized_model(description)
    if d.errors:
        return d
    reference_id = get_reference_id(d.corrected_model)
    all_transcripts = annotated_transcripts(reference_id)
    missing_transcripts = set(selectors) - all_transcripts
    if missing_transcripts:
        return {
            "errors": [{"details": f"Transcript(s) not found: {', '.join(missing_transcripts)}"}],
            "source": "conversion"
        }
    converted_descriptions = []

    if d.de_hgvs_internal_indexing_model:
        from_model = d.de_hgvs_internal_indexing_model
    else:
        from_model = None
    for selector in selectors:
        try:
            converted_model = to_hgvs_locations(
                model=from_model,
                references=d.references,
                to_coordinate_system="c",
                to_selector_id=selector,
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
    print(convert_description("NC_000011.10:g.112088971delinsTTTTT", ["NM_003002.4", "NM_014741.5"]))

    # print(overlap_transcripts(args.reference, args.start, args.end))
    # print(annotated_transcripts(args.reference))
    # print(overlap_mane_selectors(args.reference, args.start, args.end))