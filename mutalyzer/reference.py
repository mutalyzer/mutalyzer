import bisect
import copy
import re

from mutalyzer_mutator.util import reverse_complement
from mutalyzer_retriever.reference import ASSEMBLIES, get_reference_mol_type
from mutalyzer_retriever.retriever import (
    NoReferenceError,
    NoReferenceRetrieved,
    extract_feature_model,
    get_reference_model,
)

from .util import get_end, get_start, get_submodel_by_path

SELECTOR_MOL_TYPES_TYPES = ["mRNA", "ncRNA", "CDS"]
SELECTOR_FEATURE_TYPES = ["mRNA", "ncRNA", "CDS"]
COORDINATE_C_MOL_TYPES_TYPES = ["mRNA"]
COORDINATE_N_MOL_TYPES_TYPES = ["ncRNA", "transcribed RNA"]
COORDINATE_G_MOL_TYPES_TYPES = ["dna", "genomic DNA", "DNA"]
COORDINATE_P_MOL_TYPES_TYPES = ["CDS", "unassigned"]


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

    Args:
        r_m: Reference model.
        shift: Value to be subtracted from locations.
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
    if r_m["annotations"].get("qualifiers") is None:
        r_m["annotations"]["qualifiers"] = {}
    if r_m["annotations"].get("id"):
        r_m["annotations"]["qualifiers"]["chromosome_number"] = r_m["annotations"]["id"]
    if (
            r_m["annotations"].get("features")
            and len(r_m["annotations"]["features"]) >= 1
            and r_m["annotations"]["features"][0].get("qualifiers")
            and r_m["annotations"]["features"][0]["qualifiers"].get("assembly_name")
    ):
        r_m["annotations"]["qualifiers"]["assembly_name"] = r_m[
            "annotations"]["features"][0]["qualifiers"].get("assembly_name")
    _update_ensembl_ids(f_m)
    r_m["annotations"]["features"] = [f_m]
    r_m["annotations"]["id"] = f_id
    r_m["annotations"]["qualifiers"]["mol_type"] = "genomic DNA"
    location_offset = r_m["annotations"]["location"]["start"]["position"]
    update_locations(r_m["annotations"], location_offset)
    r_m["annotations"]["qualifiers"]["location_offset"] = location_offset
    return r_m


def retrieve_reference(reference_id, selector_id=None):
    try:
        r_m = get_reference_model(re.sub(r"\s+", "", reference_id), selector_id)
    except NoReferenceRetrieved:
        return None, None
    except NoReferenceError as e:
        return None, e
    if reference_id.startswith("ENS"):
        r_m = _fix_ensembl(copy.deepcopy(r_m), reference_id)
    return r_m, None


def get_feature_path(r_m, f_id, path=None):
    if path is None:
        path = []
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
    """Check if feature is on the reverse strand."""
    location = feature.get("location")
    return location and location.get("strand") == -1


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


def _is_flat_coding_gene(feature):
    """
    True for a spliceless gene->CDS gene (mitochondrial DNA, bacteria):
    one CDS child, no mRNA/ncRNA, can act as its own c. selector.
    """
    if feature.get("type") != "gene" or not feature.get("features"):
        return False
    cds_children = [f for f in feature["features"] if f.get("type") == "CDS"]
    other_children = [f for f in feature["features"] if f.get("type") in ("mRNA", "ncRNA")]
    return len(cds_children) == 1 and not other_children


def get_selector_feature(feature_model, feature_id):
    """
    Extract the feature model corresponding to the feature_id that
    can act as a selector.
    """
    for sub_feature_model in yield_feature_models(feature_model):
        if (
            sub_feature_model.get("id")
            and sub_feature_model["id"] == feature_id
            and (
                sub_feature_model.get("type") in SELECTOR_FEATURE_TYPES
                or _is_flat_coding_gene(sub_feature_model)
            )
        ):
            return sub_feature_model
    return None


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
    return None


def get_selector_feature_model(feature_model, feature_id, path=None):
    """
    Extract the feature model corresponding to the feature_id that
    can act as a selector, along with the path to reach it.
    Returns: (sub_feature, path) tuple or (None, None)
    where path is a list of indices like [0, 2, 1] to navigate features.
    """
    if path is None:
        path = []

    for i, child in enumerate(feature_model.get("features", [])):
        if (
                child.get("id")
                and child["id"] == feature_id
                and (
                    child.get("type") in SELECTOR_FEATURE_TYPES
                    or _is_flat_coding_gene(child)
                )
        ):
            return child, path + [i]

        result, result_path = get_selector_feature_model(child, feature_id, path + [i])
        if result:
            return result, result_path

    return None, None

def get_value_by_path(annotations, path):
    """
    Extract a value from the annotations given an index path for the features.

    Args:
        annotations: The annotations model structure.
        path: List of indices like [0, 1, 2] to navigate through features.

    Returns:
        The feature at the specified path, or None if not found.
    """
    current = annotations

    for i in path:
        features = current.get("features")
        if features and isinstance(features, list) and 0 <= i < len(features):
            current = features[i]
        else:
            return None

    return current


def _find_ancestor_by_type(annotations, path, feature_type):
    """
    Find the first ancestor of a given type by traversing up the path.

    Args:
        annotations: Root annotation structure
        path: List of indices
        feature_type: Type to search for (e.g., "gene", "mRNA")

    Returns:
        The ancestor feature or None
    """
    for i in range(len(path) - 1, -1, -1):
        ancestor = get_value_by_path(annotations, path[:i])
        if ancestor and ancestor.get("type") == feature_type:
            return ancestor
    return None


def _resolve_cds_and_exon_features(annotations, feature_model, path):
    if feature_model["type"] == "CDS":
        parent = get_value_by_path(annotations, path[:-1]) if path else None
        return feature_model, parent
    return _get_cds_id(feature_model), feature_model


def _set_related_feature_id(output, key, feature_model, related_feature):
    if related_feature and related_feature is not feature_model:
        output[key] = related_feature["id"]


def _set_translation_qualifiers(output, cds_feature):
    if not (cds_feature and cds_feature.get("qualifiers")):
        return
    qualifiers = cds_feature["qualifiers"]
    if qualifiers.get("translation_exception"):
        output["translation_exception"] = qualifiers["translation_exception"]
    if qualifiers.get("exception"):
        output["exception"] = qualifiers["exception"]
    if qualifiers.get("translation_table"):
        output["translation_table"] = qualifiers["translation_table"]


def get_internal_selector_model(annotations, selector_id, fix_exon=False):
    """
    Extract and flatten the selector model.

    Args:
        annotations: Root annotation structure containing genes and transcripts.
        selector_id: ID of the selector (gene, mRNA, ncRNA, or CDS).
        fix_exon: If True, creates a single exon spanning the entire transcript when none exist.

    Returns:
        dict: Flattened selector model with keys:
            - id: Selector ID.
            - type: Feature type (gene, mRNA, ncRNA, CDS).
            - gene_id: Parent (or, for a gene selector, own) gene ID.
            - exon: List of (start, end) tuples.
            - cds: List of (start, end) tuples (if coding).
            - inverted: Boolean indicating reverse strand.
        None: If selector not found.
    """
    feature_model, path = get_selector_feature_model(annotations, selector_id)

    if not feature_model:
        return None

    output = {
        "id": selector_id,
        "type": feature_model["type"],
        "inverted": is_feature_inverted(feature_model),
        "location": feature_model["location"],
    }

    # Unlike cds_id/mrna_id below a gene selector has no other,
    # different feature to point at.
    if feature_model["type"] == "gene":
        output["gene_id"] = feature_model["id"]
    else:
        gene_model = _find_ancestor_by_type(annotations, path, "gene")
        if gene_model and gene_model.get("id"):
            output["gene_id"] = gene_model["id"]

    if (
            feature_model.get("qualifiers")
            and feature_model["qualifiers"].get("tag")
            and "MANE" in feature_model["qualifiers"]["tag"]
    ):
        output["tag"] = feature_model["qualifiers"]["tag"]

    if feature_model.get("qualifiers") and feature_model["qualifiers"].get("biotype"):
        output["biotype"] = feature_model["qualifiers"]["biotype"]

    cds_feature, exon_feature = _resolve_cds_and_exon_features(
        annotations, feature_model, path
    )
    _set_related_feature_id(output, "cds_id", feature_model, cds_feature)
    _set_translation_qualifiers(output, cds_feature)
    _set_related_feature_id(output, "mrna_id", feature_model, exon_feature)
    if exon_feature:
        output.update(sort_locations(get_feature_locations(exon_feature)))

    if fix_exon and output.get("exon") is None:
        output["exon"] = [(get_start(output), get_end(output))]
        output["whole_exon_transcript"] = True

    return output


def get_available_selectors(reference_annotations, coordinate_system):
    return get_selectors_ids(reference_annotations, coordinate_system)


def get_protein_selector_model(reference, selector_id):
    selector_model = get_internal_selector_model(reference, selector_id, True)
    mrna = get_selector_feature(reference, selector_id)
    if mrna["type"] == "CDS":
        selector_model["protein_id"] = selector_id
        selector_model["transcript_id"] = selector_id
        return selector_model
    protein_ids = set()
    if mrna.get("features"):
        for feature in mrna["features"]:
            if feature["type"] == "CDS":
                protein_ids.add(feature["id"])
    if len(protein_ids) == 1:
        selector_model["protein_id"] = list(protein_ids)[0]
        selector_model["transcript_id"] = selector_id
        return selector_model
    return None


def extract_reference_id(references):
    if (
        references.get("reference")
        and references["reference"].get("model")
        and references["reference"]["model"].get("id")
    ):
        return references["reference"]["model"]["id"]
    return None

def extract_sequences(references):
    """
    Return a dictionary with reference ids as keys and their corresponding
    sequences as values.

    Args:
        references: Dictionary with reference models.
    Returns:
        Dict with Reference ids as keys and their corresponding sequences as values
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
    raise Exception("No reference ID found in the model.")


def is_selector_in_reference(selector_id, model):
    for reference_selector_id in yield_selector_ids(model):
        if selector_id == reference_selector_id:
            return True
    return False


def yield_selectors(model):
    for gene in yield_gene_models(model):
        if _is_flat_coding_gene(gene):
            yield gene
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

    Args:
        model: Reference annotations model.
        l_min: 5' location.
        l_max: 3' location.
    Returns:
        Minimum and maximum locations based on the overlapping selectors.
    """
    new_min = l_min
    new_max = l_max
    for gene in yield_gene_models(model):
        for selector in gene.get("features", []):
            if selector["type"] in SELECTOR_MOL_TYPES_TYPES:
                if is_overlap(selector, new_min, new_max):
                    if get_start(selector["location"]) < l_min:
                        new_min = get_start(selector["location"])
                    if l_max < get_end(selector["location"]):
                        new_max = get_end(selector["location"])
    return new_min, new_max


def yield_overlap_ids(model, start, end):
    for gene in yield_gene_models(model):
        if _is_flat_coding_gene(gene) and is_overlap(gene, start, end):
            yield gene
        if gene.get("features"):
            for selector in gene["features"]:
                if selector["type"] in SELECTOR_MOL_TYPES_TYPES:
                    if is_overlap(selector, start, end):
                        yield selector


def get_first_selector_id(model):
    """Get the first selector ID, or None if no selectors exist."""
    return next(yield_selector_ids(model), None)


def is_only_one_selector(model):
    """Check if the model contains exactly one selector."""
    count = 0
    for _ in yield_selectors(model):
        count += 1
        if count > 1:
            return False
    return count == 1


def _gene_selector_ids(gene):
    """
    Selector IDs a gene resolves to by name/HGNC lookup: its own id for a
    flat coding gene (mitochondrial DNA, bacteria), else its mRNA/ncRNA/CDS.
    """
    if _is_flat_coding_gene(gene):
        return [gene["id"]]
    return [
        selector["id"]
        for selector in gene.get("features", [])
        if selector["type"] in SELECTOR_MOL_TYPES_TYPES
    ]


def get_gene_selectors(gene_name, model):
    """
    Get all selector IDs for a gene identified by gene name.

    Args:
        gene_name: Gene identifier (e.g., 'SDHD').
        model: Reference model.

    Returns:
        List of selector IDs, empty list if gene not found.
    """
    for gene in yield_gene_models(model):
        if gene.get("id") == gene_name:
            return _gene_selector_ids(gene)
    return []


def get_gene_selectors_hgnc(hgnc_id, model):
    """
    Get all selector IDs for a gene identified by HGNC ID.

    Args:
        hgnc_id: HGNC identifier (e.g., '10683').
        model: Reference model.

    Returns:
        List of selector IDs, empty list if gene not found.
    """
    for gene in yield_gene_models(model):
        if gene.get("qualifiers", {}).get("HGNC") == hgnc_id:
            return _gene_selector_ids(gene)
    return []


def coordinate_system_from_mol_type(mol_type):
    if mol_type in COORDINATE_G_MOL_TYPES_TYPES:
        return "g"
    if mol_type in COORDINATE_C_MOL_TYPES_TYPES:
        return "c"
    if mol_type in COORDINATE_N_MOL_TYPES_TYPES:
        return "n"
    if mol_type in COORDINATE_P_MOL_TYPES_TYPES:
        return "p"
    return None


def get_coordinate_system_from_selector_id(model, selector_id):
    selector = get_selector_feature(model["annotations"], selector_id)
    if selector.get("type") == "gene":
        return "c"
    return coordinate_system_from_mol_type(selector.get("type"))


def get_coordinate_system_from_reference(reference):
    c_s_m = coordinate_system_from_mol_type(get_reference_mol_type(reference))

    if (
        c_s_m == "g"
        and reference["annotations"].get("qualifiers")
        and reference["annotations"]["qualifiers"].get("genome") == "mitochondrion"
    ):
        return "m"
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

    Args:
        dict model: Reference model.
        str selector_id: Id of the selector containing the slice locations.
        int strand: Reverse complement the sequence if selector is inverted.
        bool include_cds: Slice according to the CDS.
    Returns:
        Sequence slice as str.
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
    """
    All locations for the selector_id and the feature type to which they correspond.
    """
    for feature in get_selector_feature(r_model["annotations"], selector_id)[
        "features"
    ]:
        if feature.get("location"):
            yield feature["location"], feature["type"]


def get_mane_transcript(suggestions, reference_id):
    """
    Get MANE Select or MANE Plus Clinical transcript if available.
    MANE Select is preferred over MANE Plus Clinical.

    Returns:
        dict with transcript info (including 'id' and 'tag') if MANE transcript found, else None.
    """
    if not suggestions or reference_id not in suggestions:
        return None

    transcripts = suggestions[reference_id]

    mane_plus_clinical = None

    for transcript in transcripts:
        tag = transcript.get("tag")
        if tag == "MANE Select":
            return transcript
        if tag == "MANE Plus Clinical":
            mane_plus_clinical = transcript

    return mane_plus_clinical


def get_assembly_from_chr_id(chr_id):
    for a in ASSEMBLIES:
        if chr_id in ASSEMBLIES[a].values():
            return a
    return None
