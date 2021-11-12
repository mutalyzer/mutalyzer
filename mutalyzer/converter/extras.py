from Bio.SeqUtils import seq1, seq3

from ..description_model import location_to_description, yield_values
from ..util import create_exact_point_model, set_by_path
from .to_hgvs_coordinates import (
    crossmap_to_hgvs_setup,
    locations_to_hgvs_locations,
    point_to_hgvs,
)


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
