from copy import deepcopy

from Bio.SeqUtils import seq1, seq3
from mutalyzer_crossmapper import NonCoding

from ..description_model import location_to_description, yield_values
from ..util import (
    create_exact_point_model,
    get_end,
    get_inserted_sequence,
    get_start,
    set_by_path,
    set_end,
    set_start,
    slice_seq
)
from .to_hgvs_coordinates import crossmap_to_hgvs_setup, point_to_hgvs
from .to_rna import get_location_type
from ..reference import extract_feature_model, get_internal_selector_model, yield_locations_selector_id

def convert_tuples(t, x, inverted=False):
    t_c = [
        (
            location_to_description(point_to_hgvs(create_exact_point_model(e[0]), **x)),
            location_to_description(
                point_to_hgvs(create_exact_point_model(e[1] - 1), **x)
            ),
        )
        for e in t
    ]
    if inverted:
        t_c = [(e[1], e[0]) for e in t_c[::-1]]
    return t_c


def convert_selector_model(s_m):
    if s_m.get("type") == "mRNA":
        x = crossmap_to_hgvs_setup("g", s_m)
        exon_g = convert_tuples(s_m["exon"], x, s_m["inverted"])
        cds_g = convert_tuples(s_m["cds"], x, s_m["inverted"])

        x = crossmap_to_hgvs_setup("c", s_m, True)
        exon_c = convert_tuples(s_m["exon"], x, s_m["inverted"])
        cds_c = convert_tuples(s_m["cds"], x, s_m["inverted"])
        return {"exon": {"g": exon_g, "c": exon_c}, "cds": {"g": cds_g, "c": cds_c}}

    x = crossmap_to_hgvs_setup("g", s_m)
    exon_g = convert_tuples(s_m["exon"], x, s_m["inverted"])

    x = crossmap_to_hgvs_setup("n", s_m)
    exon_n = convert_tuples(s_m["exon"], x, s_m["inverted"])

    return {"exon": {"g": exon_g, "n": exon_n}}


def convert_amino_acids(model, to):
    for sequence, path in yield_values(model, ["sequence", "amino_acid"]):
        if to == "1a":
            seq_1a = str(seq1(sequence))
            if sequence == str(seq3(seq_1a)):
                set_by_path(model, path, seq_1a)
            else:
                pass
        if to == "3a":
            seq_3a = str(seq3(sequence))
            if sequence == str(seq1(seq_3a)):
                set_by_path(model, path, seq_3a)
            else:
                pass


def convert_to_exons(variants, exons, sequences):
    """
    Get variants with locations converted according to the exon locations.

    :param variants: Variant models to be converted.
    :param exons: Exons locations as list of tuples.
    :param sequences: Dictionary with sequences.
    :return: The converted and unconverted variants.
    """
    converted = []
    skipped = []
    x = NonCoding(exons).coordinate_to_noncoding
    for i, v in enumerate(variants):
        slice_v = deepcopy(v)
        if v.get("location"):
            if get_location_type(v["location"], exons, 0, 0) in [
                "same exon",
                "exon exon",
            ]:
                set_start(slice_v["location"], x(get_start(slice_v))[0] - 1)
                set_end(
                    slice_v["location"],
                    x(get_end(slice_v))[0] + x(get_end(slice_v))[1] - 1,
                )
                if slice_v.get("inserted"):
                    slice_v["inserted"] = [
                        {
                            "source": "description",
                            "sequence": get_inserted_sequence(slice_v, sequences),
                        }
                    ]
                converted.append(slice_v)
            else:
                skipped.append(i)
    return converted, skipped


def get_gene_locations(r_model):
    g_l = r_model["annotations"]["features"][0]["location"]
    return g_l["start"]["position"], g_l["end"]["position"]


def convert_reference_model(reference_model, selector_id=None, slice_to=None):
    """
    Modify the reference model based on the `selector_id`
    and the `slice_to` parameters.

    :param reference_model:
    :param selector_id:
    :param slice_to:
    :return:
    """
    if selector_id is None:
        return reference_model

    new_r_model = {
        "annotations": deepcopy(
            extract_feature_model(reference_model["annotations"], selector_id)[0]
        )
    }
    ref_seq = reference_model["sequence"]["seq"]
    s_model = get_internal_selector_model(new_r_model["annotations"], selector_id, True)

    if slice_to == "transcript":
        ref_seq = slice_seq(ref_seq, s_model["exon"])
        new_r_model["sequence"] = {"seq": ref_seq}
        x = NonCoding(s_model["exon"]).coordinate_to_noncoding
        exon_end = [exon[1] for exon in s_model["exon"]]
    if slice_to == "gene":
        g_l = get_gene_locations(new_r_model)
        ref_seq = slice_seq(ref_seq, [g_l])
        new_r_model["sequence"] = {"seq": ref_seq}
        x = NonCoding([g_l]).coordinate_to_noncoding
        exon_end = [g_l[1]]

    for loc, feature_type in yield_locations_selector_id(new_r_model, selector_id):
        loc["start"]["position"] = x(loc["start"]["position"])[0] - 1
        if feature_type == "exon" and loc["end"]["position"] in exon_end:
            loc["end"]["position"] = x(loc["end"]["position"])[0]
        else:
            loc["end"]["position"] = x(loc["end"]["position"])[0] - 1
    return new_r_model
