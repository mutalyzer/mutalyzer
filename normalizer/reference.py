import json
from functools import lru_cache
from mutalyzer_retriever import retrieve_model
from pathlib import Path


from .util import get_end, get_start, cache_dir


@lru_cache(maxsize=32)
def get_reference_model(reference_id):
    cache = cache_dir()
    if cache and (Path(cache) / reference_id).is_file():
        with open(Path(cache) / reference_id) as json_file:
            return json.load(json_file)
    return retrieve_model(reference_id)


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
    if reference_model.get("features"):
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
        output = {"id": selector_id, "type": feature["type"], "inverted": is_feature_inverted(feature)}
        output.update(sort_locations(get_feature_locations(feature)))
        return output


def get_available_selectors(reference_annotations, coordinate_system):
    return get_selectors_ids(reference_annotations, coordinate_system)


def get_protein_selector_models(reference):
    """

    :param reference: Reference annotations model (not the sequence).
    :return:
    """
    selector_ids = get_selectors_ids(reference, "c")
    for selector_id in selector_ids[:20]:
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


def point_within_feature(point, feature):
    if feature.get("location"):
        if (
            feature["location"]["start"]["position"]
            <= point
            <= feature["location"]["end"]["position"]
        ):
            return True
    return False


def get_selectors_overlap(point, reference_model):
    for feature in reference_model["features"]:
        if feature["type"] == "gene" and feature.get("features"):
            for sub_feature in feature["features"]:
                if sub_feature.get("type") and sub_feature["type"] in [
                    "mRNA",
                    "lnc_RNA",
                    "ncRNA",
                ]:
                    if point_within_feature(point, sub_feature):
                        output = {
                            "type": sub_feature["type"],
                            "id": sub_feature["id"],
                            "coordinate_system": "c"
                            if sub_feature["type"] in ["mRNA"]
                            else "n",
                            "inverted": is_feature_inverted(sub_feature),
                        }
                        output.update(
                            sort_locations(get_feature_locations(sub_feature))
                        )
                        yield output


def get_only_selector(reference_model, coordinate_system=None):
    available_selectors = get_selectors_ids(reference_model, coordinate_system)
    if len(available_selectors) == 1:
        return get_selector_model(reference_model, available_selectors[0])


def coordinate_system_from_mol_type(mol_type):
    if mol_type in ['dna', 'genomic DNA']:
        return 'g'
    elif mol_type in ['mRNA']:
        return 'c'
    elif mol_type in ['ncRNA', 'transcribed RNA']:
        return 'n'
    else:
        return ''


def coordinate_system_from_reference(reference):
    mol_type = get_mol_type(reference)
    return coordinate_system_from_mol_type(mol_type)


def coordinate_system_from_selector(selector_model):
    if selector_model["type"] in ["mRNA"]:
        return "c"
    elif selector_model["type"] in ["ncRNA"]:
        return "n"
    else:
        return ""


def get_reference_id_from_model(model):
    if model.get('annotations') and model['annotations'].get('id'):
        return model['annotations']['id']
    else:
        raise Exception('No reference ID found in the model.')


def is_selector_in_reference(selector_id, model):
    pass



class Reference(object):
    def __init__(self, reference_model):
        self.model = reference_model
        self.id = self.get_id()

    def get_selector_model(self, selector_id):
        if self.model:
            return get_selector_model(self.model["model"], selector_id)

    def get_mol_type(self):
        return get_mol_type(self.model)

    def get_only_selector(self, coordinate_system=None):
        return get_only_selector(self.model["model"], coordinate_system)

    def get_default_coordinate_system(self):
        return coordinate_system_from_reference(self.model)

    def get_available_selectors(self):
        return get_selectors_ids(self.model["model"])

    def get_length(self):
        return len(self.model["sequence"]["seq"])

    def get_id(self):
        return self.model['model']['id']

    def sequence(self):
        return self.model["sequence"]["seq"]
