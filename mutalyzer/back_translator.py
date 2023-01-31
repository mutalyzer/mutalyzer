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
    d.assembly_checks()
    d.retrieve_references()
    d.pre_conversion_checks()
    if d.corrected_model.get("type") == "description_protein":
        d.normalize_protein()

    return d.back_translated_descriptions