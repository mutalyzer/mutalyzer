import bisect
import copy
import json
import re
from functools import lru_cache
from pathlib import Path

from mutalyzer_mutator.util import reverse_complement
from mutalyzer_retriever import retrieve_model
from mutalyzer_retriever.retriever import NoReferenceError, NoReferenceRetrieved

from .util import cache_dir, get_end, get_start, get_submodel_by_path

SELECTOR_MOL_TYPES_TYPES = ["mRNA", "ncRNA", "CDS"]
SELECTOR_FEATURE_TYPES = ["mRNA", "ncRNA", "CDS"]
COORDINATE_C_MOL_TYPES_TYPES = ["mRNA"]
COORDINATE_N_MOL_TYPES_TYPES = ["ncRNA", "transcribed RNA"]
COORDINATE_G_MOL_TYPES_TYPES = ["dna", "genomic DNA", "DNA"]
COORDINATE_P_MOL_TYPES_TYPES = ["CDS", "unassigned"]


@lru_cache(maxsize=32)
def get_reference_model(r_id):
    cache = cache_dir()
    if cache and (Path(cache) / r_id).is_file():
        # print(f" - from cache {r_id}")
        with open(Path(cache) / r_id) as json_file:
            return json.load(json_file)
    # print(f" - not from cache {r_id}")
    return retrieve_model(r_id, timeout=10)


def _update_ensembl_ids(r_m):
    """
    Add the version in the id.
    """
    if r_m.get("id") and r_m.get("qualifiers") and r_m["qualifiers"].get("version"):
        r_m["id"] = r_m["id"] + "." + r_m["qualifiers"]["version"]
    if r_m.get("features"):
        for feature in r_m["features"]:
            _update_ensembl_ids(feature)


def update_locations(r_m, shift):
    """
    Update the locations of all the features in the model by subtracting
    the shift value.

    :param r_m: Reference model.
    :param shift: Value to be subtracted from locations.
    """
    if r_m.get("location"):
        if r_m["location"].get("start") and r_m["location"]["start"].get("position"):
            r_m["location"]["start"]["position"] -= shift
        if r_m["location"].get("end") and r_m["location"]["end"].get("position"):
            r_m["location"]["end"]["position"] -= shift
    if r_m.get("features"):
        for feature in r_m["features"]:
            update_locations(feature, shift)


def _fix_ensembl(r_m, r_id):
    if "." in r_id:
        r_id = r_id.split(".")[0]
    f_m = extract_feature_model(r_m["annotations"], r_id, ancestors=False)[0]
    if f_m["location"]["strand"] == -1:
        r_m["sequence"]["seq"] = reverse_complement(r_m["sequence"]["seq"])
    f_id = f_m["id"] + "." + f_m["qualifiers"]["version"]
    if f_m["type"] == "mRNA":
        f_p = get_feature_path(r_m["annotations"], r_id)
        gene_model = get_submodel_by_path(r_m["annotations"], f_p[:-2])
        gene_model["features"] = [f_m]
        f_m = gene_model
    _update_ensembl_ids(f_m)
    r_m["annotations"]["features"] = [f_m]
    r_m["annotations"]["id"] = f_id
    if r_m["annotations"].get("qualifiers") is None:
        r_m["annotations"]["qualifiers"] = {}
    r_m["annotations"]["qualifiers"]["mol_type"] = "genomic DNA"
    update_locations(
        r_m["annotations"], r_m["annotations"]["location"]["start"]["position"]
    )
    return r_m


def retrieve_reference(reference_id):
    try:
        r_m = get_reference_model(re.sub("\s+", "", reference_id))
    except NoReferenceRetrieved:
        return None, None
    except NoReferenceError as e:
        return None, e
    if reference_id.startswith("ENS"):
        r_m = _fix_ensembl(copy.deepcopy(r_m), reference_id)
    return r_m, None


def get_reference_model_segmented(
    reference_id, feature_id=None, siblings=False, ancestors=True, descendants=True
):
    reference_model = retrieve_reference(reference_id)[0]
    if feature_id is not None:
        return extract_feature_model(
            reference_model["annotations"],
            feature_id,
            siblings,
            ancestors,
            descendants,
        )[0]
    return reference_model


def extract_feature_model(
    feature, feature_id, siblings=False, ancestors=True, descendants=True
):
    output_model = None
    just_found = False
    if feature.get("id") is not None and feature_id == feature["id"]:
        output_model = copy.deepcopy(feature)
        if not descendants:
            if output_model.get("features"):
                output_model.pop("features")
            return output_model, True, True
        return output_model, True, False
    elif feature.get("features"):
        for f in feature["features"]:
            output_model, just_found, propagate = extract_feature_model(
                f, feature_id, siblings, ancestors, descendants
            )
            if output_model:
                break
        if output_model and just_found:
            if siblings:
                output_model = copy.deepcopy(feature["features"])
            if not ancestors:
                return output_model, False, True
        elif propagate:
            return output_model, False, True
    if output_model is not None:
        if isinstance(output_model, dict):
            output_model = [output_model]
        return (
            {
                **{
                    k: copy.deepcopy(feature[k])
                    for k in list(set(feature.keys()) - {"features"})
                },
                **{"features": output_model},
            },
            False,
            False,
        )
    return None, False, False


def get_feature_path(r_m, f_id, path=[]):
    r = None
    if r_m.get("id") == f_id:
        return path
    if r_m.get("features"):
        for i, f in enumerate(r_m["features"]):
            r = get_feature_path(f, f_id, path + ["features", i])
            if r:
                break
    return r


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


def get_selector_feature(feature_model, feature_id):
    """
    Extract the feature model corresponding to the feature_id that
    can act as a selector.
    """
    for sub_feature_model in yield_feature_models(feature_model):
        if (
            sub_feature_model.get("id")
            and sub_feature_model["id"] == feature_id
            and sub_feature_model.get("type") in SELECTOR_FEATURE_TYPES
        ):
            return sub_feature_model


def get_feature_locations(feature):
    sub_features_locations = {}
    if feature.get("features"):
        for sub_feature in feature["features"]:
            if sub_feature["type"].lower() not in sub_features_locations:
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


def _get_cds_id(feature_model):
    """
    Get any CDS sub feature contained in the feature model. Assumes there is
    only one such feature, otherwise it will return the first instance only.
    """
    for sub_feature_model in yield_feature_models(feature_model, False):
        if sub_feature_model.get("type") == "CDS":
            return sub_feature_model


def get_internal_selector_model(reference_annotations, selector_id, fix_exon=False):
    """
    Searches for the appropriate selector model:
    - exons and cds for coding selectors;
    - only the exons for the non-coding ones.
    The model includes the selector type.
    :return: Dictionary.
    """
    feature_model = get_selector_feature(reference_annotations, selector_id)

    if feature_model:
        output = {
            "id": selector_id,
            "type": feature_model["type"],
            "inverted": is_feature_inverted(feature_model),
            "location": feature_model["location"],
        }
        cds_sub_feature_model = _get_cds_id(feature_model)
        if cds_sub_feature_model:
            output["cds_id"] = cds_sub_feature_model["id"]
            if cds_sub_feature_model.get("qualifiers"):
                if cds_sub_feature_model["qualifiers"].get("translation_exception"):
                    output["translation_exception"] = cds_sub_feature_model[
                        "qualifiers"
                    ]["translation_exception"]
                if cds_sub_feature_model["qualifiers"].get("exception"):
                    output["exception"] = cds_sub_feature_model["qualifiers"][
                        "exception"
                    ]
        if feature_model["type"] == "CDS":
            parent_model = get_feature_parent(
                reference_annotations, feature_model["id"]
            )
            output["mrna_id"] = parent_model["id"]
            output.update(sort_locations(get_feature_locations(parent_model)))
        else:
            output.update(sort_locations(get_feature_locations(feature_model)))
        if fix_exon and output.get("exon") is None:
            output["exon"] = [(get_start(output), get_end(output))]
            output["whole_exon_transcript"] = True
        return output


def get_available_selectors(reference_annotations, coordinate_system):
    return get_selectors_ids(reference_annotations, coordinate_system)


def get_protein_selector_model(reference, selector_id):
    selector_model = get_internal_selector_model(reference, selector_id, True)
    mrna = get_selector_feature(reference, selector_id)
    protein_ids = set()
    if mrna.get("features"):
        for feature in mrna["features"]:
            if feature["type"] == "CDS":
                protein_ids.add(feature["id"])
    if len(protein_ids) == 1:
        selector_model["protein_id"] = list(protein_ids)[0]
        selector_model["transcript_id"] = selector_id
        return selector_model


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
        if selector["type"] == "mRNA" and selector.get("features"):
            for s in selector["features"]:
                if s["type"] == "CDS":
                    yield s["id"]


def yield_selector_ids_coordinate_system(model, coordinate_system):
    for selector in yield_selectors(model):
        if coordinate_system_from_mol_type(selector.get("type")) == coordinate_system:
            yield selector["id"]


def yield_feature_models(feature_model, include_top=True):
    if include_top:
        yield feature_model
    if feature_model.get("features"):
        for sub_feature in feature_model["features"]:
            yield from yield_feature_models(sub_feature)


def get_feature_parent(feature_model, feature_id):
    path = get_feature_path(feature_model, feature_id)
    return get_submodel_by_path(feature_model, path[:-2])


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


def get_only_selector_id(model):
    for selector_id in yield_selector_ids(model):
        return selector_id


def is_only_one_selector(model):
    i = 0
    for selector_id in yield_selectors(model):
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
    elif mol_type in COORDINATE_P_MOL_TYPES_TYPES:
        return "p"
    return None


def get_coordinate_system_from_selector_id(model, selector_id):
    selector = get_selector_feature(model["annotations"], selector_id)
    return coordinate_system_from_mol_type(selector.get("type"))


def get_reference_mol_type(model):
    return get_model_qualifier(model, "mol_type")


def get_coordinate_system_from_reference(reference):
    c_s_m = coordinate_system_from_mol_type(get_reference_mol_type(reference))

    if (
        c_s_m == "g"
        and reference["annotations"].get("qualifiers")
        and reference["annotations"]["qualifiers"].get("genome") == "mitochondrion"
    ):
        return "m"
    else:
        return c_s_m


def _get_exons_and_cds(s_m):
    exons = [e for t in s_m["exon"] for e in t]
    cds = [s_m["cds"][0][0], s_m["cds"][-1][1]]
    return exons, cds


def _get_cds_into_exons(s_m):
    exons, cds = _get_exons_and_cds(s_m)
    l_index = bisect.bisect_right(exons, cds[0])
    r_index = bisect.bisect_left(exons, cds[1])
    new_exons = [cds[0]] + exons[l_index:r_index] + [cds[1]]
    return list(zip(new_exons[0::2], new_exons[1::2]))


def slice_to_selector(model, selector_id, strand=False, include_cds=False):
    """
    Slice the reference model sequence according to the exons and cds
    locations of the selector with the provided id.

    :arg dict model: Reference model.
    :arg str selector_id: Id of the selector containing the slice locations.
    :arg int strand: Reverse complement the sequence if selector is inverted.
    :arg bool include_cds: Slice according to the CDS.
    :returns: Sequence slice.
    :rtype: str
    """
    s_m = get_internal_selector_model(model["annotations"], selector_id, True)
    output = ""
    slices = s_m["exon"]
    if include_cds and s_m.get("cds"):
        slices = _get_cds_into_exons(s_m)
    for s in slices:
        output += model["sequence"]["seq"][s[0] : s[1]]
    if strand and s_m["inverted"]:
        output = reverse_complement(output)
    return output


def yield_locations(annotations):
    """
    All locations present in the annotations and the feature to which they
    correspond.
    """
    if annotations.get("location"):
        yield annotations["location"], annotations["type"]
    if annotations.get("features"):
        for feature in annotations["features"]:
            yield from yield_locations(feature)


def yield_locations_selector_id(r_model, selector_id):
    for feature in get_selector_feature(r_model["annotations"], selector_id)[
        "features"
    ]:
        if feature.get("location"):
            yield feature["location"], feature["type"]


GRCH38 = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "23": "NC_000023.11",
    "X": "NC_000023.11",
    "24": "NC_000024.10",
    "Y": "NC_000024.10",
}

GRCH37 = {
    "1": "NC_000001.10",
    "2": "NC_000002.11",
    "3": "NC_000003.11",
    "4": "NC_000004.11",
    "5": "NC_000005.9",
    "6": "NC_000006.11",
    "7": "NC_000007.13",
    "8": "NC_000008.10",
    "9": "NC_000009.11",
    "10": "NC_000010.10",
    "11": "NC_000011.9",
    "12": "NC_000012.11",
    "13": "NC_000013.10",
    "14": "NC_000014.8",
    "15": "NC_000015.9",
    "16": "NC_000016.9",
    "17": "NC_000017.10",
    "18": "NC_000018.9",
    "19": "NC_000019.9",
    "20": "NC_000020.10",
    "21": "NC_000021.8",
    "22": "NC_000022.10",
    "23": "NC_000023.10",
    "X": "NC_000023.10",
    "24": "NC_000024.9",
    "Y": "NC_000024.9",
}

ASSEMBLY_ALIASES = {
    "HG38": "GRCH38",
    "HG19": "GRCH37",
}

ASSEMBLIES = {
    "GRCH38": GRCH38,
    "GRCH37": GRCH37,
}


def get_model_qualifier(model, qualifier):
    if (
        model.get("annotations")
        and model["annotations"].get("qualifiers")
        and model["annotations"]["qualifiers"]
    ):
        return model["annotations"]["qualifiers"].get(qualifier)


def get_chromosome_accession(ref_id, model):
    if get_reference_mol_type(model) == "mRNA" and ref_id.startswith("NM"):
        chromosome_number = get_model_qualifier(model, "chromosome")
        if chromosome_number is not None:
            chromosome_accessions = []
            for assembly in ASSEMBLIES:
                chromosome_accession = ASSEMBLIES[assembly].get(chromosome_number)
                if chromosome_accession is not None:
                    chromosome_accessions.append((assembly, chromosome_accession))
            if chromosome_accessions:
                return chromosome_accessions
