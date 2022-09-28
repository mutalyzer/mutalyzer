from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from mutalyzer_crossmapper import Coding
from mutalyzer_mutator import mutate
from mutalyzer_mutator.util import reverse_complement

from .converter import to_rna_protein_coordinates
from .reference import extract_sequences
from .util import get_start


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
        return "p.(=)", 0, 0, 0

    lcp = len(longest_common_prefix(s1, s2))
    lcs = len(longest_common_suffix(s1[lcp:], s2[lcp:]))
    s1_end = len(s1) - lcs
    s2_end = len(s2) - lcs

    # Insertion / Duplication / Extension.
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


def get_protein_sequence(reference_model, selector_model):
    exons = selector_model["exon"]
    cds = [selector_model["cds"][0][0], selector_model["cds"][0][1]]
    dna_ref_seq = reference_model["sequence"]["seq"]
    cds_seq = slice_seq(dna_ref_seq, exons, cds[0], cds[1])
    if selector_model["inverted"]:
        cds_seq = reverse_complement(cds_seq)
    seq = list(str(Seq(cds_seq).translate()))
    if selector_model.get("translation_exception"):
        x = Coding(
            selector_model["exon"], selector_model["cds"][0], selector_model["inverted"]
        )
        for t_e in selector_model.get("translation_exception")["exceptions"]:
            seq[x.coordinate_to_protein(get_start(t_e))[0] - 1] = t_e["amino_acid"]
    return "".join(seq)


def get_protein_references(references, selector_model):
    pass


def get_protein_description(variants, references, selector_model):
    """
    Retrieves the protein description.

    :param variants: Only deletion_insertion variants with coordinate locations.
                     Preferable, from the description extractor.
    :param references: References models. Required to be able to retrieve the
                       inserted sequences.
    :param selector_model: The selector model that includes the exon and cds
                           information.
    """
    sequences = extract_sequences(references)
    ref_id = references["reference"]["annotations"]["id"]
    dna_ref_seq = sequences[ref_id]
    exons = selector_model["exon"]
    cds = [selector_model["cds"][0][0], selector_model["cds"][0][1]]
    protein_id = selector_model["protein_id"]

    cds_seq = slice_seq(dna_ref_seq, exons, cds[0], cds[1])

    if selector_model["inverted"]:
        cds_seq = reverse_complement(cds_seq)
        cds_seq_ext = reverse_complement(slice_seq(dna_ref_seq, exons, 0, cds[1]))
    else:
        cds_seq_ext = slice_seq(dna_ref_seq, exons, cds[0])

    p_ref_seq = str(Seq(cds_seq).translate())

    cds_variants, splice_site_hits = to_rna_protein_coordinates(
        variants, sequences, selector_model
    )

    if splice_site_hits:
        return "{}({}):{}".format(ref_id, protein_id, "p.?"), p_ref_seq, "?"
    elif not cds_variants:
        return "{}({}):{}".format(ref_id, protein_id, "p.(=)"), p_ref_seq, p_ref_seq

    cds_obs_seq = mutate({"reference": cds_seq_ext}, cds_variants)

    p_obs_seq = str(Seq(cds_obs_seq).translate())

    if cds_seq[:3] != cds_obs_seq[:3]:
        return "{}({}):{}".format(ref_id, protein_id, "p.?"), p_ref_seq, "?"

    # Up to and including the first '*', or the entire string.
    try:
        stop = p_obs_seq.index("*")
        p_obs_seq = p_obs_seq[: stop + 1]
    except ValueError:
        pass

    cds_stop = len(mutate({"reference": cds_seq}, cds_variants))
    description = protein_description(cds_stop, p_ref_seq, p_obs_seq)
    if len(cds_variants) > 1 and "*" in description[0]:
        # TODO: This seems to happen in M2. Check why.
        # A different check maybe should be implemented
        # see: NG_012337.1(NM_012459.2):c.5_6delinsTAG
        return (
            "{}({}):{}".format(ref_id, protein_id, "p.?"),
            p_ref_seq,
            p_obs_seq,
        )

    return (
        "{}({}):{}".format(ref_id, protein_id, description[0]),
        p_ref_seq,
        p_obs_seq,
    )
