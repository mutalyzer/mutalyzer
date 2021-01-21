import json

from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from mutalyzer_mutator import mutate
from mutalyzer_mutator.util import reverse_complement

from .converter import to_cds_coordinate
from .reference import (
    extract_sequences,
    get_protein_selector_models,
    get_reference_mol_type,
)


def longest_common_prefix(s1, s2):
    """
    Calculate the longest common prefix of two strings.

        >>> longest_common_prefix('abcdefg', 'abcabcdefg')
        'abc'
        >>> longest_common_prefix('abcdefg', 'abcdefg')
        'abcdefg'

    @arg s1: The first string.
    @type s1: unicode
    @arg s2: The second string.
    @type s2: unicode

    @return: The longest common prefix of s1 and s2.
    @rtype: unicode

    @todo: This is mostly used just for the length of the returned string,
           and we could also return that directly.
    """
    pos = 0

    while pos < min(len(s1), len(s2)) and s1[pos] == s2[pos]:
        pos += 1

    return s1[:pos]


def longest_common_suffix(s1, s2):
    """
    Calculate the longest common suffix of two strings.

        >>> longest_common_suffix('abcdefg', 'abcabcdefg')
        'abcdefg'
        >>> longest_common_suffix('abcdefg', 'abcefg')
        'efg'

    @arg s1: The first string.
    @type s1: unicode
    @arg s2: The second string.
    @type s2: unicode

    @return: The longest common suffix of s1 and s2.
    @rtype: unicode
    """
    return longest_common_prefix(s1[::-1], s2[::-1])[::-1]


def in_frame_description(s1, s2):
    """
    Give a description of an inframe difference of two proteins. Also give
    the position at which the proteins start to differ and the positions at
    which they are the same again.

        >>> in_frame_description('MTAPQQMT*', 'MTAQQMT*')
        ('p.(Pro4del)', 3, 4, 3)
        >>> in_frame_description('MTAPQQMT*', 'MTAQMT*')
        ('p.(Pro4_Gln5del)', 3, 5, 3)
        >>> in_frame_description('MTAPQQT*', 'MTAQQMT*')
        ('p.(Pro4_Gln6delinsGlnGlnMet)', 3, 6, 6)
        >>> in_frame_description('MTAPQQMT*', 'MTAPQQMTMQ*')
        ('p.(*9Metext*2)', 8, 9, 11)
        >>> in_frame_description('MTAPQQMT*', 'MTAPQQMTMQ')
        ('p.(*9Metext*?)', 8, 9, 10)

    @arg s1: The original protein.
    @type s1: unicode
    @arg s2: The mutated protein.
    @type s2: unicode

    @return: A tuple of:
        - unicode ; Protein description of the change.
        - int     ; First position of the change.
        - int     ; Last position of the change in the first protein.
        - int     ; Last position of the change in the second protein.
    @rtype: tuple(unicode, int, int, int)

    @todo: More intelligently handle longest_common_prefix().
    @todo: Refactor this code (too many return statements).
    """
    s2_stop = "*" in s2
    s1 = s1.rstrip("*")
    s2 = s2.rstrip("*")

    if s1 == s2:
        # Nothing happened.
        return ("p.(=)", 0, 0, 0)

    lcp = len(longest_common_prefix(s1, s2))
    lcs = len(longest_common_suffix(s1[lcp:], s2[lcp:]))
    s1_end = len(s1) - lcs
    s2_end = len(s2) - lcs

    # Insertion / Duplication / Extention.
    if not s1_end - lcp:
        if len(s1) == lcp:
            # http://www.hgvs.org/mutnomen/FAQ.html#nostop
            stop = str(abs(len(s1) - len(s2))) if s2_stop else "?"

            return (
                "p.(*%i%sext*%s)" % (len(s1) + 1, seq3(s2[len(s1)]), stop),
                len(s1),
                len(s1) + 1,
                len(s2) + (1 if s2_stop else 0),
            )

        ins_length = s2_end - lcp

        if lcp - ins_length >= 0 and s1[lcp - ins_length : lcp] == s2[lcp:s2_end]:
            if ins_length == 1:
                return (
                    "p.(%s%idup)" % (seq3(s1[lcp - ins_length]), lcp - ins_length + 1),
                    lcp,
                    lcp,
                    lcp + 1,
                )
            return (
                "p.(%s%i_%s%idup)"
                % (
                    seq3(s1[lcp - ins_length]),
                    lcp - ins_length + 1,
                    seq3(s1[lcp - 1]),
                    lcp,
                ),
                lcp,
                lcp,
                lcp + ins_length,
            )
        # if
        return (
            "p.(%s%i_%s%iins%s)"
            % (seq3(s1[lcp - 1]), lcp, seq3(s1[lcp]), lcp + 1, seq3(s2[lcp:s2_end])),
            lcp,
            lcp,
            s2_end,
        )
    # if

    # Deletion / Inframe stop.
    if not s2_end - lcp:
        if len(s2) == lcp:
            return (
                "p.(%s%i*)" % (seq3(s1[len(s2)]), len(s2) + 1),
                lcp,
                len(s1) + 1,
                len(s2) + 1,
            )

        if lcp + 1 == s1_end:
            return ("p.(%s%idel)" % (seq3(s1[lcp]), lcp + 1), lcp, lcp + 1, lcp)
        return (
            "p.(%s%i_%s%idel)" % (seq3(s1[lcp]), lcp + 1, seq3(s1[s1_end - 1]), s1_end),
            lcp,
            s1_end,
            lcp,
        )
    # if

    # Substitution.
    if s1_end == s2_end and s1_end == lcp + 1:
        return (
            "p.(%s%i%s)" % (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp])),
            lcp,
            lcp + 1,
            lcp + 1,
        )

    # InDel.
    if lcp + 1 == s1_end:
        return (
            "p.(%s%idelins%s)" % (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp:s2_end])),
            lcp,
            lcp + 1,
            s2_end,
        )
    return (
        "p.(%s%i_%s%idelins%s)"
        % (seq3(s1[lcp]), lcp + 1, seq3(s1[s1_end - 1]), s1_end, seq3(s2[lcp:s2_end])),
        lcp,
        s1_end,
        s2_end,
    )


def out_of_frame_description(s1, s2):
    """
    Give the description of an out of frame difference between two
    proteins. Give a description of an inframe difference of two proteins.
    Also give the position at which the proteins start to differ and the
    end positions (to be compatible with the in_frame_description function).

        >>> out_of_frame_description('MTAPQQMT*', 'MTAQQMT*')
        ('p.(Pro4Glnfs*5)', 3, 9, 8)
        >>> out_of_frame_description('MTAPQQMT*', 'MTAQMT*')
        ('p.(Pro4Glnfs*4)', 3, 9, 7)
        >>> out_of_frame_description('MTAPQQT*', 'MTAQQMT*')
        ('p.(Pro4Glnfs*5)', 3, 8, 8)
        >>> out_of_frame_description('MTAPQQT*', 'MTAQQMT')
        ('p.(Pro4Glnfs*?)', 3, 8, 7)

    @arg s1: The original protein.
    @type s1: unicode
    @arg s2: The mutated protein.
    @type s2: unicode

    @return: A tuple of:
        - unicode ; Protein description of the change.
        - int     ; First position of the change.
        - int     ; Last position of the first protein.
        - int     ; Last position of the second protein.
    @rtype: tuple(unicode, int, int, int)

    @todo: More intelligently handle longest_common_prefix().
    """
    s1_seq = s1.rstrip("*")
    s2_seq = s2.rstrip("*")
    lcp = len(longest_common_prefix(s1_seq, s2_seq))

    if lcp == len(s2_seq):  # NonSense mutation.
        if lcp == len(s1_seq):  # Is this correct?
            return ("p.(=)", 0, 0, 0)
        return ("p.(%s%i*)" % (seq3(s1[lcp]), lcp + 1), lcp, len(s1), lcp)
    if lcp == len(s1_seq):
        # http://www.hgvs.org/mutnomen/FAQ.html#nostop
        stop = str(abs(len(s1_seq) - len(s2_seq))) if "*" in s2 else "?"

        return (
            "p.(*%i%sext*%s)" % (len(s1_seq) + 1, seq3(s2[len(s1_seq)]), stop),
            len(s1_seq),
            len(s1),
            len(s2),
        )

    # http://www.hgvs.org/mutnomen/FAQ.html#nostop
    stop = str(len(s2_seq) - lcp + 1) if "*" in s2 else "?"

    return (
        "p.(%s%i%sfs*%s)" % (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp]), stop),
        lcp,
        len(s1),
        len(s2),
    )


def protein_description(cds_stop, s1, s2):
    """
    Wrapper function for the in_frame_description() and
    out_of_frame_description() functions. It uses the value cds_stop to
    decide which one to call.

        >>> protein_description(34, 'MTAPQQMT*', 'MTAQQMT*')
        ('p.(Pro4Glnfs*5)', 3, 9, 8)
        >>> protein_description(34, 'MTAPQQMT*', 'MTAQQMT')
        ('p.(Pro4Glnfs*?)', 3, 9, 7)
        >>> protein_description(33, 'MTAPQQMT*', 'MTAQQMT*')
        ('p.(Pro4del)', 3, 4, 3)
        >>> protein_description(33, 'MTAPQQMT*', 'TTAQQMT*')
        ('p.?', 0, 4, 3)

    @arg cds_stop: Position of the stop codon in c. notation (CDS length).
    @type cds_stop: int
    @arg s1: The original protein.
    @type s1: unicode
    @arg s2: The mutated protein.
    @type s2: unicode

    @return: A tuple of:
        - unicode ; Protein description of the change.
        - int     ; First position of the change.
        - int     ; Last position of the change in the first protein.
        - int     ; Last position of the change in the second protein.
    @rtype: tuple(unicode, int, int, int)
    """
    if cds_stop % 3:
        description = out_of_frame_description(s1, s2)
    else:
        description = in_frame_description(s1, s2)

    return description


def new_index(value, slices):
    output = 0
    for s in slices:
        if s[0] <= value <= s[1]:
            output += value - s[0]
            break
        else:
            output += s[1] - s[0]
    return output


def slice_seq(seq, slices, start=None, end=None):
    output = ""
    for s in slices:
        output += seq[s[0] : s[1]]
    start = new_index(start, slices) if start else 0
    end = new_index(end, slices) if end else -1
    return output[start:end]


def get_protein_description(variants, references, selector_model):
    """
    Retrieves the protein description.

    :param variants: Only deletion_insertion variants with coordinate locations.
                     Preferable, from the description extractor.
    :param references: References models. Required to be able to retrieve the
                       inserted sequences.
    :param selector_model:
    """
    sequences = extract_sequences(references)

    cds_variants = to_cds_coordinate(variants, sequences, selector_model)

    cds_sequence = slice_seq(
        sequences[references["reference"]["annotations"]["id"]],
        selector_model["exon"],
        selector_model["cds"][0][0],
        selector_model["cds"][0][1],
    )

    if selector_model["inverted"]:
        cds_sequence = reverse_complement(cds_sequence)
        cds_sequence_extended = reverse_complement(
            slice_seq(
                sequences[references["reference"]["annotations"]["id"]],
                selector_model["exon"],
                0,
                selector_model["cds"][0][1],
            )
        )
    else:
        cds_sequence_extended = slice_seq(
            sequences[references["reference"]["annotations"]["id"]],
            selector_model["exon"],
            selector_model["cds"][0][0],
        )

    cds_sequence_mutated = mutate({"reference": cds_sequence_extended}, cds_variants)

    reference_protein = Seq(cds_sequence).translate()
    predicted_protein = Seq(cds_sequence_mutated).translate()

    if cds_sequence[:3] != cds_sequence_mutated[:3]:
        return "{}({}):{}".format(
            references["reference"]["annotations"]["id"],
            selector_model["protein_id"],
            "p.?",
        )

    # Up to and including the first '*', or the entire string.
    try:
        stop = str(predicted_protein).index("*")
        predicted_protein = str(predicted_protein)[: stop + 1]
    except ValueError:
        pass

    cds_stop = len(mutate({"reference": cds_sequence}, cds_variants))
    description = protein_description(
        cds_stop, str(reference_protein), str(predicted_protein)
    )

    return "{}({}):{}".format(
        references["reference"]["annotations"]["id"],
        selector_model["protein_id"],
        description[0],
    )


def get_protein_descriptions(variants, references):
    """
    Retrieves all the possible protein descriptions.

    - First, it identifies if protein descriptions can be obtained for this
      particular reference.
    - Next, it retrieves the selector models for which a protein description
      can be obtained.
    - Finally, it loops over the selector models and calls the appropriate
      function to obtain the protein descriptions.

    :param variants: Only deletion_insertion variants with coordinate locations.
                     Preferable, from the description extractor.
    :param references: References models. Required to be able to retrieve the
                       inserted sequences.
    """
    if get_reference_mol_type(references["reference"]) not in [
        "genomic DNA",
        "mRNA",
        "dna",
    ]:
        return

    protein_descriptions = []
    for selector_model in get_protein_selector_models(
        references["reference"]["annotations"]
    ):
        protein_descriptions.append(
            get_protein_description(variants, references, selector_model)
        )
    return protein_descriptions
