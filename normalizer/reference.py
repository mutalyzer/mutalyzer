import json

from .util import get_end, get_start


def get_mol_type(reference):
    if reference["model"].get("qualifiers"):
        if reference["model"]["qualifiers"].get("mol_type"):
            return reference["model"]["qualifiers"]["mol_type"]


def is_feature_inverted(feature):
    if feature.get("location") and feature["location"].get("strand"):
        if feature["location"]["strand"] == -1:
            return True
        else:
            return False
    else:
        # TODO: to be fixed by checking the reference model.
        return False


def get_selectors_ids(reference_annotations, coordinate_system=None):
    ids = set()
    if coordinate_system is None:
        check = ["mRNA", "ncRNA"]
    elif coordinate_system == "c":
        check = ["mRNA"]
    elif coordinate_system == "n":
        check = ["ncRNA"]
    else:
        check = []
        # TODO: raise some error.
    if reference_annotations.get("features"):
        for feature in reference_annotations["features"]:
            if feature["type"] == "gene" and feature.get("features"):
                for sub_feature in feature["features"]:
                    if sub_feature["type"] in check:
                        ids.add(sub_feature["id"])
    return list(ids)


def is_id_equal(feature, feature_id):
    """
    Runs a series of checks to identify if the feature has the provided ID.
    """
    if feature_id == feature["id"]:
        return True
    # if '-' in feature['id'] and feature_id == feature['id'].split('-')[1]:
    #     return True
    return False


def get_feature(reference_model, feature_id):
    """
    Extract the feature model, if found, otherwise None.
    """
    for feature in reference_model["features"]:
        if feature["type"] == "gene" and feature.get("features"):
            for sub_feature in feature["features"]:
                if is_id_equal(sub_feature, feature_id):
                    return sub_feature


def get_feature_locations(feature):
    sub_features_locations = {}
    if feature.get("features"):
        for sub_feature in feature["features"]:
            if sub_feature["type"] not in sub_features_locations:
                sub_features_locations[sub_feature["type"].lower()] = []
            sub_features_locations[sub_feature["type"].lower()].append(
                (get_start(sub_feature), get_end(sub_feature))
            )
    return sub_features_locations


def sort_locations(locations):
    sorted_locations = {}
    for k in locations:
        sorted_locations[k] = sorted(locations[k], key=lambda i: i[1])
    return sorted_locations


def get_selector_model(reference_model, selector_id):
    """
    Searches for the appropriate selector model:
    - exons and cds for coding selectors;
    - only the exons for the non-coding ones.
    The model includes the selector type.
    :return: Dictionary.
    """
    feature = get_feature(reference_model, selector_id)
    if feature:
        output = {"type": feature["type"], "inverted": is_feature_inverted(feature)}
        output.update(sort_locations(get_feature_locations(feature)))
        return output


def get_available_selectors(reference_annotations, coordinate_system):
    return get_selectors_ids(reference_annotations, coordinate_system)


def get_protein_selector_models(reference):
    """

    :param reference: Reference annotations model (not the sequence).
    :return:
    """
    selector_models = {}
    selector_ids = get_selectors_ids(reference, "c")
    for selector_id in selector_ids:
        selector_model = get_selector_model(reference, selector_id)
        mrna = get_feature(reference, selector_id)
        protein_ids = []
        if mrna.get("features"):
            for feature in mrna["features"]:
                if feature["type"] == "CDS":
                    protein_ids.append(feature["id"])
        if len(protein_ids) == 1:
            selector_model["protein_id"] = protein_ids[0]
            selector_model["transcript_id"] = selector_id
            selector_models[selector_id] = selector_model

    return selector_models


def extract_reference_id(references):
    if (
        references.get("reference")
        and references["reference"].get("model")
        and references["reference"]["model"].get("id")
    ):
        return references["reference"]["model"]["id"]


def extract_sequences(references):
    """
    Return a dictionary with reference ids as keys and their corresponding
    sequences as values.

    :param references: Dictionary with reference models.
    :rtype: dict
    :return: Reference ids as keys and their corresponding sequences as values
    """
    sequences = {}
    for reference in references:
        sequences[reference] = references[reference]["sequence"]["seq"]
    return sequences


def get_sequence_length(references, reference_id):
    return len(references[reference_id]["sequence"]["seq"])
