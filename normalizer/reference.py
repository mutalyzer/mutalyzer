import copy
import json
from functools import lru_cache
from pathlib import Path

from mutalyzer_retriever import retrieve_model

from .util import cache_dir, get_end, get_start

SELECTOR_MOL_TYPES_TYPES = ["mRNA", "ncRNA"]
COORDINATE_C_MOL_TYPES_TYPES = ["mRNA"]
COORDINATE_N_MOL_TYPES_TYPES = ["ncRNA", "transcribed RNA"]
COORDINATE_G_MOL_TYPES_TYPES = ["dna", "genomic DNA", "DNA"]


@lru_cache(maxsize=32)
def get_reference_model(reference_id):
    cache = cache_dir()
    if cache and (Path(cache) / reference_id).is_file():
        with open(Path(cache) / reference_id) as json_file:
            return json.load(json_file)
    return retrieve_model(reference_id, timeout=10)


def get_reference_model_segmented(
    reference_id, feature_id=None, include_siblings=False
):
    reference_model = get_reference_model(reference_id)
    if feature_id is not None:
        return extract_feature_model(
            reference_model["annotations"], feature_id, include_siblings
        )[0]
    return reference_model


def extract_feature_model(feature, feature_id, include_siblings=False):
    output_model = None
    just_found = False
    if feature.get("id") is not None and feature_id == feature["id"]:
        return copy.deepcopy(feature), True
    elif feature.get("features"):
        for f in feature["features"]:
            out = extract_feature_model(f, feature_id, include_siblings)
            output_model, just_found = out
            if output_model:
                break
        if output_model and just_found and include_siblings:
            output_model = copy.deepcopy(feature["features"])
    if output_model is not None:
        if isinstance(output_model, dict):
            output_model = [output_model]
        return {
            **{
                k: copy.deepcopy(feature[k])
                for k in list(set(feature.keys()) - {"features"})
            },
            **{"features": output_model},
        }, False
    return None, False


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


def get_feature(reference_annotations, feature_id):
    """
    Extract the feature model, if found, otherwise None.
    """
    if reference_annotations.get("features"):
        for feature in reference_annotations["features"]:
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


def get_selector_model(reference_annotations, selector_id, fix_exon=False):
    """
    Searches for the appropriate selector model:
    - exons and cds for coding selectors;
    - only the exons for the non-coding ones.
    The model includes the selector type.
    :return: Dictionary.
    """
    feature = get_feature(reference_annotations, selector_id)
    if feature:
        output = {
            "id": selector_id,
            "type": feature["type"],
            "inverted": is_feature_inverted(feature),
            "location": feature["location"],
        }
        output.update(sort_locations(get_feature_locations(feature)))
        if fix_exon and output.get("exon") is None:
            output["exon"] = [(get_start(output), get_end(output))]
        return output


def get_available_selectors(reference_annotations, coordinate_system):
    return get_selectors_ids(reference_annotations, coordinate_system)


def get_protein_selector_model(reference, selector_id):
    selector_model = get_selector_model(reference, selector_id, True)
    mrna = get_feature(reference, selector_id)
    protein_ids = []
    if mrna.get("features"):
        for feature in mrna["features"]:
            if feature["type"] == "CDS":
                protein_ids.append(feature["id"])
    if len(protein_ids) == 1:
        selector_model["protein_id"] = protein_ids[0]
        selector_model["transcript_id"] = selector_id
        return selector_model


def get_protein_selector_models(reference):
    """

    :param reference: Reference annotations model (not the sequence).
    :return:
    """
    selector_ids = get_selectors_ids(reference, "c")
    for selector_id in selector_ids[:20]:
        selector_model = get_selector_model(reference, selector_id, True)
        mrna = get_feature(reference, selector_id)
        protein_ids = []
        if mrna.get("features"):
            for feature in mrna["features"]:
                if feature["type"] == "CDS":
                    protein_ids.append(feature["id"])
        if len(protein_ids) == 1:
            selector_model["protein_id"] = protein_ids[0]
            selector_model["transcript_id"] = selector_id
            yield selector_model


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


def get_reference_id_from_model(model):
    if model.get("annotations") and model["annotations"].get("id"):
        return model["annotations"]["id"]
    else:
        raise Exception("No reference ID found in the model.")


def is_selector_in_reference(selector_id, model):
    for reference_selector_id in yield_selector_ids(model):
        if selector_id == reference_selector_id:
            return True
    return False


def yield_selectors(model):
    for gene in yield_gene_models(model):
        if gene.get("features"):
            for selector in gene["features"]:
                if selector["type"] in SELECTOR_MOL_TYPES_TYPES:
                    yield selector


def yield_gene_models(model):
    annotations = model["annotations"]
    if annotations.get("features"):
        for feature in annotations["features"]:
            if feature["type"] == "gene":
                yield feature


def yield_selector_ids(model):
    for selector in yield_selectors(model):
        yield selector["id"]


def yield_selector_ids_coordinate_system(model, coordinate_system):
    for selector in yield_selectors(model):
        if coordinate_system_from_mol_type(selector.get("type")) == coordinate_system:
            yield selector["id"]


def is_overlap(selector, start, end):
    sel_s = get_start(selector["location"])
    sel_e = get_end(selector["location"])
    if (
        sel_s <= start <= sel_e
        or sel_s <= end <= sel_e
        or start <= sel_s <= sel_e <= end
    ):
        return True
    return False


def overlap_min_max(model, l_min, l_max):
    """
    Get the overlapping minimum and maximum locations based on the selectors
    that are contain the l_min and l_max locations.

    :param model: Reference annotations model.
    :param l_min: 5' location.
    :param l_max: 3' location.
    :return: Minimum and maximum locations based on the overlapping selectors.
    """
    new_min = l_min
    new_max = l_max
    for gene in yield_gene_models(model):
        if gene.get("features"):
            for selector in gene["features"]:
                if selector["type"] in SELECTOR_MOL_TYPES_TYPES:
                    if is_overlap(selector, new_min, new_max):
                        if get_start(selector["location"]) < l_min:
                            new_min = get_start(selector["location"])
                        if l_max < get_end(selector["location"]):
                            new_max = get_end(selector["location"])
    return new_min, new_max


def yield_overlap_ids(model, start, end):
    for gene in yield_gene_models(model):
        if gene.get("features"):
            for selector in gene["features"]:
                if selector["type"] in SELECTOR_MOL_TYPES_TYPES:
                    if is_overlap(selector, start, end):
                        yield selector


def get_only_selector_id(model, coordinate_system):
    for selector_id in yield_selector_ids_coordinate_system(model, coordinate_system):
        return selector_id


def is_only_one_selector(model, coordinate_system):
    i = 0
    for selector_id in yield_selector_ids_coordinate_system(model, coordinate_system):
        if i > 1:
            break
        i += 1
    if i == 1:
        return True
    else:
        return False


def get_gene_selectors(gene_name, model):
    selectors = []
    for gene in yield_gene_models(model):
        if gene.get("id") == gene_name and gene.get("features"):
            for selector in gene["features"]:
                if selector["type"] in SELECTOR_MOL_TYPES_TYPES:
                    selectors.append(selector["id"])
            break
    return selectors


def get_gene_selectors_hgnc(hgnc_id, model):
    selectors = []
    for gene in yield_gene_models(model):
        if (
            gene.get("qualifiers")
            and gene["qualifiers"].get("HGNC") == hgnc_id
            and gene.get("features")
        ):
            for selector in gene["features"]:
                if selector["type"] in SELECTOR_MOL_TYPES_TYPES:
                    selectors.append(selector["id"])
            break
    return selectors


def coordinate_system_from_mol_type(mol_type):
    if mol_type in COORDINATE_G_MOL_TYPES_TYPES:
        return "g"
    elif mol_type in COORDINATE_C_MOL_TYPES_TYPES:
        return "c"
    elif mol_type in COORDINATE_N_MOL_TYPES_TYPES:
        return "n"
    return None


def get_coordinate_system_from_selector_id(model, selector_id):
    selector = get_feature(model["annotations"], selector_id)
    return coordinate_system_from_mol_type(selector.get("type"))


def get_reference_mol_type(model):
    if model["annotations"].get("qualifiers"):
        if model["annotations"]["qualifiers"].get("mol_type"):
            return model["annotations"]["qualifiers"]["mol_type"]


def get_coordinate_system_from_reference(reference):
    mol_type = get_reference_mol_type(reference)
    return coordinate_system_from_mol_type(mol_type)
