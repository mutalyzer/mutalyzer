"""
Back translation (p. to c.).
"""

import itertools

from mutalyzer_backtranslate import BackTranslate
from mutalyzer_mutator.util import reverse_complement

from .description import Description
from .description_model import get_reference_id
from .protein import slice_seq
from .util import get_end, get_start


def back_translate(description):
    """
    Back translation of an amino acid substitution to all possible
    causal one-nucleotide substitutions.

    :arg unicode description: Amino acid substitution in HGVS format.
    :returns: List of DNA descriptions in HGVS format.
    :rtype: list
    """
    d = Description(description)
    d.retrieve_references()
    d.pre_conversion_checks()
    # TODO: check if protein description

    bt = BackTranslate()
    translated_vars = []

    if d.get_selector_model():
        selector_model = d.get_selector_model()
        exons = selector_model["exon"]
        cds = [selector_model["cds"][0][0], selector_model["cds"][0][1]]
        dna_ref_seq = d.references["reference"]["sequence"]["seq"]
        cds_seq = slice_seq(dna_ref_seq, exons, cds[0], cds[1])

        d.normalize()

        if selector_model["inverted"]:
            cds_seq = reverse_complement(cds_seq)

        for variant in d.internal_indexing_model["variants"]:
            if variant.get("type") == "substitution":
                cds_start = get_start(variant) * 3
                cds_end = get_end(variant) * 3
                bt_options = bt.with_dna(
                    cds_seq[cds_start:cds_end],
                    variant["inserted"][0]["sequence"],
                )
                dna_var_options = []
                for offset in bt_options:
                    for v in bt_options[offset]:
                        dna_var_options.append(
                            "{}{}>{}".format(cds_start + offset + 1, v[0], v[1])
                        )
                translated_vars.append(dna_var_options)
            else:
                # TODO: Add error message.
                return []
        if selector_model.get("id") == get_reference_id(d.corrected_model):
            reference = "{}".format(get_reference_id(d.corrected_model))
        else:
            reference = "{}({})".format(
                get_reference_id(d.corrected_model), selector_model["mrna_id"]
            )
    else:
        # TODO: Should be reached for "NP_" references: get the "NM_".
        return []
    translated_descriptions = []
    for t in itertools.product(*translated_vars):
        if len(t) > 1:
            translated_descriptions.append("{}:c.([{}])".format(reference, ";".join(t)))
        elif len(t) == 1:
            translated_descriptions.append("{}:c.({})".format(reference, t[0]))

    return translated_descriptions
