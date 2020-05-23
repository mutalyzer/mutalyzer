from Bio.SeqUtils import seq3


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
#longest_common_prefix


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
#longest_common_suffix


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
    s2_stop = '*' in s2
    s1 = s1.rstrip('*')
    s2 = s2.rstrip('*')

    if s1 == s2:
        # Nothing happened.
        return ('p.(=)', 0, 0, 0)

    lcp = len(longest_common_prefix(s1, s2))
    lcs = len(longest_common_suffix(s1[lcp:], s2[lcp:]))
    s1_end = len(s1) - lcs
    s2_end = len(s2) - lcs

    # Insertion / Duplication / Extention.
    if not s1_end - lcp:
        if len(s1) == lcp:
            # http://www.hgvs.org/mutnomen/FAQ.html#nostop
            stop = str(abs(len(s1) - len(s2))) if s2_stop else '?'

            return ('p.(*%i%sext*%s)' % \
                    (len(s1) + 1, seq3(s2[len(s1)]), stop),
                    len(s1), len(s1) + 1, len(s2) + (1 if s2_stop else 0))

        ins_length = s2_end - lcp

        if lcp - ins_length >= 0 and s1[lcp - ins_length:lcp] == s2[lcp:s2_end]:
            if ins_length == 1:
                return ('p.(%s%idup)' % \
                        (seq3(s1[lcp - ins_length]), lcp - ins_length + 1),
                        lcp, lcp, lcp + 1)
            return ('p.(%s%i_%s%idup)' % \
                    (seq3(s1[lcp - ins_length]),
                     lcp - ins_length + 1, seq3(s1[lcp - 1]), lcp),
                    lcp, lcp, lcp + ins_length)
        #if
        return ('p.(%s%i_%s%iins%s)' % \
                (seq3(s1[lcp - 1]), lcp, seq3(s1[lcp]),
                 lcp + 1, seq3(s2[lcp:s2_end])),
                lcp, lcp, s2_end)
    #if

    # Deletion / Inframe stop.
    if not s2_end - lcp:
        if len(s2) == lcp:
            return ('p.(%s%i*)' % (seq3(s1[len(s2)]), len(s2) + 1),
                    lcp, len(s1) + 1, len(s2) + 1)

        if lcp + 1 == s1_end:
            return ('p.(%s%idel)' % (seq3(s1[lcp]), lcp + 1),
                    lcp, lcp + 1, lcp)
        return ('p.(%s%i_%s%idel)' % \
                (seq3(s1[lcp]), lcp + 1, seq3(s1[s1_end - 1]), s1_end),
                lcp, s1_end, lcp)
    #if

    # Substitution.
    if s1_end == s2_end and s1_end == lcp + 1:
        return ('p.(%s%i%s)' % (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp])),
                lcp, lcp + 1, lcp + 1)

    # InDel.
    if lcp + 1 == s1_end:
        return ('p.(%s%idelins%s)' % \
                (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp:s2_end])),
                lcp, lcp + 1, s2_end)
    return ('p.(%s%i_%s%idelins%s)' % \
            (seq3(s1[lcp]), lcp + 1, seq3(s1[s1_end - 1]), s1_end,
             seq3(s2[lcp:s2_end])),
            lcp, s1_end, s2_end)
#in_frame_description


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
    s1_seq = s1.rstrip('*')
    s2_seq = s2.rstrip('*')
    lcp = len(longest_common_prefix(s1_seq, s2_seq))

    if lcp == len(s2_seq): # NonSense mutation.
        if lcp == len(s1_seq): # Is this correct?
            return ('p.(=)', 0, 0, 0)
        return ('p.(%s%i*)' % (seq3(s1[lcp]), lcp + 1), lcp, len(s1), lcp)
    if lcp == len(s1_seq):
        # http://www.hgvs.org/mutnomen/FAQ.html#nostop
        stop = str(abs(len(s1_seq) - len(s2_seq))) if '*' in s2 else '?'

        return ('p.(*%i%sext*%s)' % \
                (len(s1_seq) + 1, seq3(s2[len(s1_seq)]), stop),
                len(s1_seq), len(s1), len(s2))

    # http://www.hgvs.org/mutnomen/FAQ.html#nostop
    stop = str(len(s2_seq) - lcp + 1) if '*' in s2 else '?'

    return ('p.(%s%i%sfs*%s)' % \
            (seq3(s1[lcp]), lcp + 1, seq3(s2[lcp]), stop),
            lcp, len(s1), len(s2))
#out_of_frame_description


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
#protein_description




if __name__ == '__main__':
    print(protein_description(34, 'MTAPQQMT*', 'MTAQQMT*'))
