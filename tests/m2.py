TESTS = [
    {
        "keywords": [
            "M2: test_deletion_in_frame",
            "Simple in-frame deletion should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_163del",
        "input": "NG_007485.1(NM_000077.4):c.161_163del",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_163del",
                # "AL449423.14(CDKN2A_i001):p.(Met54_Gly55delinsSer)",
                "NG_007485.1(NP_000068.1):p.(Met54_Gly55delinsSer)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_206del",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGlu)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_insertion_in_frame",
            "Simple in-frame insertion should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162insATC",
        "input": "NG_007485.1(NM_000077.4):c.161_162insATC",
        "normalized": "NG_007485.1(NM_000077.4):c.161_162insATC",
        "genomic": "NG_007485.1:g.28294_28295insATC",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162insATC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsIleSer)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsIleSer)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205insATC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69insIle)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444insGAT"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_insertion_list_in_frame",
            "Simple in-frame insertion should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
            "input not for test_description_to_model_to_description",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162insATC",
        "input": "NG_007485.1(NM_000077.4):c.161_162ins[ATC]",
        "normalized": "NG_007485.1(NM_000077.4):c.161_162insATC",
        "genomic": "NG_007485.1:g.28294_28295insATC",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162insATC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsIleSer)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsIleSer)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205insATC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69insIle)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444insGAT"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_insertion_in_frame",
            "Simple in-frame deletion/insertion should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delinsATCCC",
        "input": "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
        "genomic": "NG_007485.1:g.28294_28295delinsATCCC",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444delinsGGGAT"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_insertion_list_in_frame",
            "Simple in-frame deletion-insertion of a list should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
            "input not for test_description_to_model_to_description",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delins[ATCCC]",
        "input": "NG_007485.1(NM_000077.4):c.161_162delins[ATCCC]",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444delinsGGGAT"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_insertion_in_frame_complete",
            "Simple in-frame deletion-insertion of a list should give a simple description on protein level.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delTGinsATCCC",
        "input": "NG_007485.1(NM_000077.4):c.161_162delTGinsATCCC",
        "normalized": "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
        "genomic": "NG_007485.1:g.28294_28295delinsATCCC",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444delinsGGGAT"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_insertion_list_in_frame_complete",
            "Simple in-frame deletion-insertion of a list should give a simple description on protein level, also with the optional deleted sequence argument.",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
            "input not for test_description_to_model_to_description",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.161_162delTGinsATCCC",
        "input": "NG_007485.1(NM_000077.4):c.161_162delTGins[ATCCC]",
        "normalized": "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
        "genomic": "NG_007485.1:g.28294_28295delinsATCCC",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.161_162delinsATCCC",
                # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
                "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
            ),
            (
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444delinsGGGAT"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_est_warning_nm_est",
            "M2: Warning for EST positioning on NM reference.",
            "Switch to info on coordinate system correction.",
        ],
        "input": "NM_003002.2:274del",
        "normalized": "NM_003002.2:c.274del",
        "infos": ["ICORRECTEDCOORDINATESYSTEM"],
        "to_test": True,
    },
    # test_no_est_warning_nm_c: no longer relevant ?
    # test_no_est_warning_nm_n: no longer relevant ?
    {
        "keywords": [
            "M2: test_est_warning_ng_est",
            "M2: Warning for EST positioning on NG reference.",
            "Switch to info on coordinate system correction and other reference ID.",
        ],
        "input": "NG_012337.1:128del",
        "normalized": "NG_012337.1:g.128del",
        "infos": ["ICORRECTEDCOORDINATESYSTEM"],
        "to_test": True,
    },
    # test_no_est_warning_ng_g: no longer relevant ?
    # test_no_est_warning_est_est: no longer relevant ?
    {
        "keywords": [
            "M2: test_roll",
            "M2: Just a variant where we should roll.",
        ],
        "input": "NM_003002.2:c.273del",
        "normalized": "NM_003002.2:c.274del",
        "coding_protein_descriptions": {
            ("NM_003002.2:c.274del", "NM_003002.2(NP_002993.1):p.(Asp92Thrfs*43)")
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_no_roll",
            "M2: Just a variant where we cannot roll.",
        ],
        "input": "NM_003002.2:c.274del",
        "normalized": "NM_003002.2:c.274del",
        "coding_protein_descriptions": {
            ("NM_003002.2:c.274del", "NM_003002.2(NP_002993.1):p.(Asp92Thrfs*43)")
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_no_roll_splice",
            "M2: Here we can roll but should not, because it is over a splice site.",
            # TODO: In M3 a coordinate system mismatch error is raised.
            # TODO: In any case in M3 we do roll. Should we?
        ],
        "input": "NM_000088.3:g.459del",
        "normalized": "",
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_partial_roll_splice",
            "M2: Here we can roll two positions, but should roll only one because otherwise it is over a splice site.",
            # TODO: Check previous.
        ],
        "input": "NM_000088.3:g.494del",
        "normalized": "",
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_roll_after_splice",
            "M2: Here we can roll and should, we stay in the same exon.",
            # TODO: Check previous.
        ],
        "input": "NM_000088.3:g.460del",
        "normalized": "",
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_roll_both_ins",
            """M2:
            Insertion that rolls should not use the same inserted sequence in
            descriptions on forward and reverse strands.
            Here we have the following situation on the forward strand:
                                      65470 (genomic)
                                        |
                CGGTGCGTTGGGCAGCGCCCCCGCCTCCAGCAGCGCCCGCACCTCCTCTA
            Now, an insertion of TAC after 65470 should be rolled to an insertion
            of ACT after 65471:

                CGGTGCGTTGGGCAGCGCCCCCGCC --- TCCAGCAGCGCCCGCACCTCCTCTA
                CGGTGCGTTGGGCAGCGCCCCCGCC TAC TCCAGCAGCGCCCGCACCTCCTCTA  =>

                CGGTGCGTTGGGCAGCGCCCCCGCCT --- CCAGCAGCGCCCGCACCTCCTCTA
                CGGTGCGTTGGGCAGCGCCCCCGCCT ACT CCAGCAGCGCCCGCACCTCCTCTA
            However, in CDKN2A_v001 (on the reverse strand), this insertion should
            roll the other direction and the inserted sequence should be the reverse
            complement of CTA, which is TAG, and not that of ACT, which is AGT.
            The next test (test_roll_reverse_ins) tests the situation for an input
            of AL449423.14:g.65471_65472insACT, where only the reverse roll should
            be done.
            """,
            "Switched from AL449423.14 to NG_007485.1 and updated locations.",
        ],
        # "input": "AL449423.14:g.65470_65471insTAC",
        "input": "NG_007485.1:g.5478_5479insTAC",
        "normalized": "NG_007485.1:g.5479_5480insACT",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.193+126_193+127insACT",
                "NG_007485.1(NP_478102.2):p.(=)",
            ),
            (
                "NG_007485.1(NM_000077.4):c.-19186_-19185insACT",
                "NG_007485.1(NP_000068.1):p.(=)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+26260_616+26261insTAG"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_roll_reverse_ins",
            """M2:
            Insertion that rolls on the reverse strand should not use the same
            inserted sequence in descriptions on forward and reverse strands.
            """,
            "Switched from AL449423.14 to NG_007485.1 and updated locations.",
        ],
        # "input": "AL449423.14:g.65471_65472insACT",
        "input": "NG_007485.1:g.5479_5480insACT",
        "normalized": "NG_007485.1:g.5479_5480insACT",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.193+126_193+127insACT",
                "NG_007485.1(NP_478102.2):p.(=)",
            ),
            (
                "NG_007485.1(NM_000077.4):c.-19186_-19185insACT",
                "NG_007485.1(NP_000068.1):p.(=)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+26260_616+26261insTAG"],
        "to_test": False,
    },
    # test_roll_message_forward: switch to check the shift amount?
    # test_roll_message_reverse: switch to check the shift amount?
    {
        "keywords": [
            "M2: test_ins_cds_start",
            "M2: Insertion on CDS start boundary should not be included in CDS.",
        ],
        "input": "NM_000143.3:c.-1_1insCAT",
        "normalized": "NM_000143.3:c.-1_1insCAT",
        "coding_protein_descriptions": {
            (
                "NM_000143.3:c.-1_1insCAT",
                "NM_000143.3(NP_000134.2):p.(=)",  # M2
                # TODO: or "NM_000143.3(NP_000134.2):p.?",
            ),
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_ins_cds_start_after",
            "M2: Insertion after CDS start boundary should be included in CDS.",
        ],
        "input": "NM_000143.3:c.1_2insCAT",
        "normalized": "NM_000143.3:c.1_2insCAT",
        "coding_protein_descriptions": {
            (
                "NM_000143.3:c.1_2insCAT",
                "NM_000143.3(NP_000134.2):p.?",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_splice_site",
            "M2: Deletion hitting one splice site should not do a protein prediction.",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.632-5_670del",
        "normalized": "NG_012772.1(NM_000059.3):c.632-5_670del",
        "genomic": "NG_012772.1:g.18959_19002del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.632-5_670del",
                "NG_012772.1(NP_000050.2):p.?",
            ),
        },
        # TODO: Add splice site warning?
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_del_exon",
            "M2: Deletion of an entire exon should be possible.",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.632-5_681+7del",
        "normalized": "NG_012772.1(NM_000059.3):c.632-5_681+7del",
        "genomic": "NG_012772.1:g.18959_19020del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.632-5_681+7del",
                "NG_012772.1(NP_000050.2):p.(Val211Glufs*10)",
            ),
        },
        # TODO: Add splice site warning?
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_exon_exact",
            "M2: Deletion of exactly an exon should be possible.",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.632_681del",
        "normalized": "NG_012772.1(NM_000059.3):c.632_681del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.632_681del",
                "NG_012772.1(NP_000050.2):p.(Val211Glufs*10)",
            ),
        },
        # TODO: Add splice site warning?
        "to_test": True,
    },
    {
        "keywords": [
            "M2:  test_del_exon_in_frame",
            """M2:
            Deletion of an entire exon with length a triplicate should give a
            proteine product with just this deletion (and possibly substitutions
            directly before and after).
            NG_012772.1(BRCA2_v001):c.68-7_316+7del is such a variant, since
            positions 68 through 316 are exactly one exon and (316-68+1)/3 = 83.
            """,
        ],
        "input": "NG_012772.1(BRCA2_v001):c.68-7_316+7del",
        "normalized": "NG_012772.1(NM_000059.3):c.68-7_316+7del",
        "genomic": "NG_012772.1:g.8591_8853del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.68-7_316+7del",
                "NG_012772.1(NP_000050.2):p.(Asp23_Leu105del)",
            ),
        },
        # TODO: Add splice site warning?
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_exons",
            "M2: Deletion of two entire exons should be possible.",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.632-5_793+7del",
        "normalized": "NG_012772.1(NM_000059.3):c.632-5_793+7del",
        "genomic": "NG_012772.1:g.18959_20558del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.632-5_793+7del",
                "NG_012772.1(NP_000050.2):p.(Val211_His264del)",
            ),
        },
        # TODO: Add splice site warning?
        "to_test": True,
    },
    {
        "keywords": [
            "M2:  test_del_intron",
            "M2: Deletion of an entire intron should be possible (fusion of remaining exonic parts).",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.622_674del",
        "normalized": "NG_012772.1(NM_000059.3):c.622_674del",
        "genomic": "NG_012772.1:g.16125_19006del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.622_674del",
                "NG_012772.1(NP_000050.2):p.(Val208Tyrfs*3)",
            ),
        },
        # TODO: Add splice site warning?
        "to_test": True,
    },
    {
        "keywords": [
            "M2:  test_del_intron_exact",
            "M2: Deletion of exactly an intron should be possible (fusion of flanking exons).",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.681+1_682-1del",
        "normalized": "NG_012772.1(NM_000059.3):c.681+1_682-1del",
        "genomic": "NG_012772.1:g.19014_20439del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.681+1_682-1del",
                "NG_012772.1(NP_000050.2):p.(=)",
            ),
        },
        # TODO: Add splice site warning?
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_intron_in_frame",
            "M2: Deletion of an entire intron should be possible (fusion of remaining exonic parts).",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.622_672del",
        "normalized": "NG_012772.1(NM_000059.3):c.622_672del",
        "genomic": "NG_012772.1:g.16125_19004del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.622_672del",
                "NG_012772.1(NP_000050.2):p.(Val208_Asp224del)",
            ),
        },
        # TODO: Add splice site warning?
        "to_test": True,
    },
    {
        "keywords": [
            "M2: del_exon_unknown_offsets",
            "Deletion of an entire exon with unknown offsets should be possible.",
            "To be implemented.",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.632-?_681+?del",
        "normalized": "NG_012772.1(NM_000059.3):c.632-?_681+?del",
        "genomic": "NG_012772.1:g.(17550_19725)del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.632-?_681+?del",
                "NG_012772.1(NP_000050.2):p.(Val211Glufs*10)",
            )
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_del_exon_unknown_offsets_in_frame",
            """M2:
            Deletion of an entire exon with unknown offsets and length a
            triplicate should give a proteine product with just this deletion
            (and possibly substitutions directly before and after).
            NG_012772.1(BRCA2_v001):c.68-?_316+?del is such a variant, since
            positions 68 through 316 are exactly one exon and (316-68+1)/3 = 83.
            """,
        ],
        # M2: Genomic positions should be centered in flanking introns and unsure],
        # M2: Todo: .c notation should still be c.632-?_681+?del, but what about other transcripts?
        "input": "NG_012772.1(BRCA2_v001):c.68-?_316+?del",
        "normalized": "NG_012772.1(NM_000059.3):c.68-?_316+?del",
        "genomic": "NG_012772.1:g.(7324_11720)del",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.68-?_316+?del",
                "NG_012772.1(NP_000050.2):p.(Asp23_Leu105del)",
            )
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_del_exon_unknown_offsets_composed",
            "M2: Deletion of an entire exon with unknown offsets and another "
            "composed variant with exact positioning should be possible.",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.[632-?_681+?del;681+4del]",
        "normalized": "NG_012772.1(NM_000059.3):c.[632-?_681+?del;681+4del]",
        "genomic": "NG_012772.1:g.[(17550_19725)del;19017del]",
        "coding_protein_descriptions": {
            (
                "NG_012772.1(NM_000059.3):c.[632-?_681+?del;681+4del]",
                "NG_012772.1(NP_000050.2):p.?",
            )
        },
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_del_exon_unknown_offsets_reverse",
            "M2: Deletion of an entire exon with unknown offsets should be "
            "possible, also on the reverse strand.",
            # TODO
            "To be adapted and implemented",
        ],
        "input": "AL449423.14(CDKN2A_v001):c.151-?_457+?del",
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_del_exon_transcript_reference",
            """M2:
            Deletion of entire exon on a transcript reference should remove the
            expected splice sites (only that of the deleted exon), and not those
            of the flanking exons (as would happen using the mechanism for genomic
            references).
            """,
        ],
        "input": "NM_000143.3:c.739_904del",
        "normalized": "NM_000143.3:c.739_904del",
        "coding_protein_descriptions": {
            (
                "NM_000143.3:c.739_904del",
                "NM_000143.3(NP_000134.2):p.(Glu247Alafs*27)",
            )
        },
        # TODO: Add splice site warning?
        # TODO: M2 message: Sequence "GAATTT [154bp] TTACAG" at position
        #  802_967 was not corrected to "AATTTA [154bp] TACAGG" at position
        #  803_968, since they reside in different exons.
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_ins_seq",
            "M2: Insertion of a sequence.",
        ],
        "input": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_reverse",
            "M2: Insertion of a sequence on reverse strand.",
            "Transcripts mismatch between genbank and gff3.",
        ],
        "input": "NG_012337.1(TIMM8B_v001):c.12_13insGATC",
        "normalized": "NG_012337.1(NM_012459.2):c.12_13insGATC",
        "genomic": "NG_012337.1:g.4911_4912insATCG",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_012459.2):c.12_13insGATC",
                "NG_012337.1(NP_036591.2):p.(Ser5Aspfs*21)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range",
            "M2: Insertion of a range.",
        ],
        "input": "NG_008939.1:g.5207_5208ins4300_4320",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_inv",
            "M2: Insertion of an inverse range.",
        ],
        "input": "NG_008939.1:g.5207_5208ins4300_4320inv",
        "normalized": "NG_008939.1:g.5207_5208insGCCAGATAATGAGCACAGGAC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGCCAGATAATGAGCACAGGAC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Leu539delinsAlaArg)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_list",
            "M2: Insertion of a sequence as a list.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1:g.5207_5208ins[GTCCTGTGCTCATTATCTGGC]",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_list_reverse",
            "M2: Insertion of a sequence as a list on reverse strand.",
            "Transcripts mismatch between genbank and gff3.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_012337.1(TIMM8B_v001):c.12_13ins[GATC]",
        "normalized": "NG_012337.1(NM_012459.2):c.12_13insGATC",
        "genomic": "NG_012337.1:g.4911_4912insATCG",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_012459.2):c.12_13insGATC",
                "NG_012337.1(NP_036591.2):p.(Ser5Aspfs*21)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_list",
            "M2: Insertion of a range as a list.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1:g.5207_5208ins[4300_4320]",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_inv_list",
            "M2: Insertion of an inverse range as a list.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1:g.5207_5208ins[4300_4320inv]",
        "normalized": "NG_008939.1:g.5207_5208insGCCAGATAATGAGCACAGGAC",
        "genomic": "NG_008939.1:g.5207_5208insGCCAGATAATGAGCACAGGAC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGCCAGATAATGAGCACAGGAC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Leu539delinsAlaArg)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_seq",
            "M2: Insertion of two sequences.",
        ],
        "input": "NG_008939.1:g.5207_5208ins[GTCCTGTGCTC;ATTATCTGGC]",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_seq_reverse",
            "M2: Insertion of two sequences on reverse strand.",
        ],
        "input": "NG_012337.1(TIMM8B_v001):c.12_13ins[TTT;GATC]",
        "normalized": "NG_012337.1(NM_012459.2):c.12_13insTTTGATC",
        "genomic": "NG_012337.1:g.4911_4912insATCAAAG",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_012459.2):c.12_13insTTTGATC",
                "NG_012337.1(NP_036591.2):p.(Ser5Phefs*22)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_range",
            "M2: Insertion of two ranges.",
        ],
        "input": "NG_008939.1:g.5207_5208ins[4300_4309;4310_4320]",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_range_inv",
            "M2: Insertion of a range and an inverse range.",
        ],
        "input": "NG_008939.1:g.5207_5208ins[4300_4309;4310_4320inv]",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTGCCAGATAATG",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTGCCAGATAATG",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTGCCAGATAATG",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysCysGlnIleMet)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_range",
            "M2: Insertion of a sequence and a range.",
        ],
        "input": "NG_008939.1:g.5207_5208ins[GTCCTGTGCT;4310_4320]",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_range_inv",
            "M2: Insertion of a sequence and an inverse range.",
        ],
        "input": "NG_008939.1:g.5207_5208ins[GTCCTGTGCT;4310_4320inv]",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTGCCAGATAATG",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTGCCAGATAATG",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTGCCAGATAATG",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysCysGlnIleMet)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_seq",
            "M2: Insertion of a range and a sequence.",
        ],
        "input": "NG_008939.1:g.5207_5208ins[4300_4309;CATTATCTGGC]",
        "normalized": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_inv_seq",
            "M2: Insertion of an inverse range and a sequence.",
        ],
        "input": "NG_008939.1:g.5207_5208ins[4300_4309inv;CATTATCTGGC]",
        "normalized": "NG_008939.1:g.5207_5208insAGCACAGGACCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5208insAGCACAGGACCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insAGCACAGGACCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insSerThrGlyProLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_coding",
            "M2: Insertion of a sequence (coding).",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_157insGTCCTGTGCTCATTATCTGGC",
        "normalized": "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_list_coding",
            "M2: Insertion of a sequence as a list (coding).",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_157ins[GTCCTGTGCTCATTATCTGGC]",
        "normalized": "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_seq_seq_coding",
            "M2: Insertion of two sequences (coding).",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_157ins[GTCCTGTGCTC;ATTATCTGGC]",
        "normalized": "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5208insGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157insGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_coding",
            "M2: Insertion of a range (coding).",
            "Updated, since M3 supports this insertion.",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_157ins180_188",
        "normalized": "NG_008939.1(NM_000532.5):c.156_157ins180_188",
        "genomic": "NG_008939.1:g.5207_5208ins5231_10536",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157ins180_188",
                "NG_008939.1(NP_000523.2):p.(Arg53Alafs*5)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_inv_coding",
            "M2: Insertion of an inverse range (coding).",
            "Updated, since M3 supports this insertion.",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_157ins180_188inv",
        "normalized": "NG_008939.1(NM_000532.5):c.156_157ins180_188inv",
        "genomic": "NG_008939.1:g.5207_5208ins5231_10536inv",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157ins180_188inv",
                "NG_008939.1(NP_000523.2):p.(Arg53Phefs*11)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_list_coding",
            "M2: Insertion of a range as a list (coding).",
            "Updated, since M3 supports this insertion.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_157ins[180_188]",
        "normalized": "NG_008939.1(NM_000532.5):c.156_157ins180_188",
        "genomic": "NG_008939.1:g.5207_5208ins5231_10536",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157ins180_188",
                "NG_008939.1(NP_000523.2):p.(Arg53Alafs*5)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_range_inv_list_coding",
            "M2: Insertion of an inverse range as a list (coding).",
            "Updated, since M3 supports this insertion.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_157ins[180_188inv]",
        "normalized": "NG_008939.1(NM_000532.5):c.156_157ins180_188inv",
        "genomic": "NG_008939.1:g.5207_5208ins5231_10536inv",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_157ins180_188inv",
                "NG_008939.1(NP_000523.2):p.(Arg53Phefs*11)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_seq",
            "M2: Insertion-deletion of a sequence.",
        ],
        "input": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "normalized": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range",
            "M2: Insertion-deletion of a range.",
        ],
        "input": "NG_008939.1:g.5207_5212delins4300_4320",
        "normalized": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_inv",
            "M2: Insertion-deletion of an inverse range.",
        ],
        "input": "NG_008939.1:g.5207_5212delins4300_4320inv",
        "normalized": "NG_008939.1:g.5207_5212delinsGCCAGATAATGAGCACAGGAC",
        "genomic": "NG_008939.1:g.5207_5212delinsGCCAGATAATGAGCACAGGAC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGCCAGATAATGAGCACAGGAC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsProAspAsnGluHisArgThr)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_seq_list",
            "M2: Insertion-deletion of a sequence as a list.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1:g.5207_5212delins[GTCCTGTGCTCATTATCTGGC]",
        "normalized": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_list",
            "M2: Insertion-deletion of a range as a list.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1:g.5207_5212delins[4300_4320]",
        "normalized": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_inv_list",
            "M2: Insertion-deletion of an inverse range as a list.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1:g.5207_5212delins[4300_4320inv]",
        "normalized": "NG_008939.1:g.5207_5212delinsGCCAGATAATGAGCACAGGAC",
        "genomic": "NG_008939.1:g.5207_5212delinsGCCAGATAATGAGCACAGGAC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGCCAGATAATGAGCACAGGAC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsProAspAsnGluHisArgThr)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_seq_seq",
            "M2: Insertion-deletion of two sequences.",
        ],
        "input": "NG_008939.1:g.5207_5212delins[GTCCTGTGCT;CATTATCTGGC]",
        "normalized": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_range",
            "M2: Insertion-deletion of two ranges.",
        ],
        "input": "NG_008939.1:g.5207_5212delins[4300_4309;4310_4320]",
        "normalized": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_inv_range",
            """M2:
            Insertion-deletion of an inverse range and a range.
            Note that the delins is also shortened by one position here.
            """,
        ],
        "input": "NG_008939.1:g.5207_5212delins[4300_4309inv;4310_4320]",
        "normalized": "NG_008939.1:g.5208_5212delinsGCACAGGACCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.157_161delinsGCACAGGACCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsAlaGlnAspHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_seq_range",
            "M2: Insertion-deletion of a sequence and a range.",
        ],
        "input": "NG_008939.1:g.5207_5212delins[GTCCTGTGCT;4310_4320]",
        "normalized": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_seq_range_inv",
            """M2:
            Insertion-deletion of a sequence and an inverse range.
            Note that the delins is also shortened by one position here.
            """,
        ],
        "input": "NG_008939.1:g.5207_5212delins[GTCCTGTGCT;4310_4320inv]",
        "normalized": "NG_008939.1:g.5207_5211delinsGTCCTGTGCTGCCAGATAAT",
        "genomic": "NG_008939.1:g.5207_5211delinsGTCCTGTGCTGCCAGATAAT",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_160delinsGTCCTGTGCTGCCAGATAAT",
                "NG_008939.1(NP_000523.2):p.(Arg53_Leu539delinsSerCysAlaAlaArg)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_seq",
            "M2: Insertion-deletion of a range and a sequence.",
        ],
        "input": "NG_008939.1:g.5207_5212delins[4300_4309;CATTATCTGGC]",
        "normalized": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_inv_seq",
            """M2:
            Insertion-deletion of an inverse range and a sequence.
            Note that the delins is also shortened by one position here.
            """,
        ],
        "input": "NG_008939.1:g.5207_5212delins[4300_4309inv;CATTATCTGGC]",
        "normalized": "NG_008939.1:g.5208_5212delinsGCACAGGACCATTATCTGGC",
        "genomic": "NG_008939.1:g.5208_5212delinsGCACAGGACCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.157_161delinsGCACAGGACCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsAlaGlnAspHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_seq_coding",
            "M2: Insertion-deletion of a sequence (coding).",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
        "normalized": "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_seq_list_coding",
            "M2: Insertion-deletion of a sequence as a list (coding).",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_161delins[GTCCTGTGCTCATTATCTGGC]",
        "normalized": "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_seq_seq_coding",
            "M2: Insertion-deletion of two sequences (coding).",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_161delins[GTCCTGTGCT;CATTATCTGGC]",
        "normalized": "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
        "genomic": "NG_008939.1:g.5207_5212delinsGTCCTGTGCTCATTATCTGGC",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.156_161delinsGTCCTGTGCTCATTATCTGGC",
                "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_coding",
            "M2: Insertion-deletion of a range (coding).",
            "Updated, since M3 supports this insertion.",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_161delins180_188",
        "normalized": "NG_008939.1(NM_000532.5):c.[155_156ins180_183+69;156_161inv;161_162ins183+76_188]",
        "genomic": "NG_008939.1:g.[5206_5207ins5231_5303;5207_5212inv;5212_5213ins5310_10536]",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_inv_coding",
            "M2: Insertion-deletion of an inverse range (coding).",
            "Updated, since M3 supports this insertion.",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_161delins180_188inv",
        "normalized": "NG_008939.1(NM_000532.5):c.[155_156ins183+76_188inv;161_162ins180_183+69inv]",
        "genomic": "NG_008939.1:g.[5206_5207ins5310_10536inv;5212_5213ins5231_5303inv]",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_list_coding",
            "M2: Insertion-deletion of a range as a list (coding).",
            "Updated, since M3 supports this insertion.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_161delins[180_188]",
        "normalized": "NG_008939.1(NM_000532.5):c.[155_156ins180_183+69;156_161inv;161_162ins183+76_188]",
        "genomic": "NG_008939.1:g.[5206_5207ins5231_5303;5207_5212inv;5212_5213ins5310_10536]",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_inv_list_coding",
            "M2: Insertion-deletion of an inverse range as a list (coding).",
            "Updated, since M3 supports this insertion.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "NG_008939.1(PCCB_v001):c.156_161delins[180_188inv]",
        "normalized": "NG_008939.1(NM_000532.5):c.[155_156ins183+76_188inv;161_162ins180_183+69inv]",
        "genomic": "NG_008939.1:g.[5206_5207ins5310_10536inv;5212_5213ins5231_5303inv]",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_no_reference",
            "M2: Variant description without a reference.",
            "M2 error ENOREF: No reference sequence given.",
            "input not for test_description_to_model_to_description",
        ],
        "input": "g.244355733del",
        "errors": ["ESYNTAXUEOF"],
        "to_test": True,
    },
    # test_chromosomal_positions: seems like a position convert test?
    # test_ex_notation: not supported in M3
    {
        "keywords": [
            "M2: test_lrg_reference",
            "M2: We should be able to use LRG reference sequence without error.",
        ],
        "input": "LRG_1t1:c.266G>T",
        "normalized": "LRG_1(t1):c.266G>T",  # TODO: check if OK.
        "genomic": "LRG_1:g.6855G>T",
        "coding_protein_descriptions": {
            (
                "LRG_1(t1):c.266G>T",
                "LRG_1(p1):p.(Gly89Val)",
            )
        },
        "infos": ["ICORRECTEDLRGREFERENCE"],
        "to_test": True,
    },
    # test_gi_reference_plain: gi numbers not supported in M3.
    # test_gi_reference_prefix:
    # test_gi_reference_prefix_colon:
    {
        "keywords": [
            "M2: test_nop_nm",
            "M2: Variant on NM without effect should be described as '='.",
        ],
        "input": "NM_002001.2:c.1_3delinsATG",
        "normalized": "NM_002001.2:c.=",
        "to_test": True,
    },
    # test_nop_ud: no UDs in M3.
    # test_ud_reverse_sequence:
    # test_ud_forward_sequence:
    # test_ud_reverse_range:
    # test_ud_forward_range:
    # test_ud_reverse_del_length:
    # test_ud_reverse_roll:
    # test_ud_forward_roll:
    {
        "keywords": [
            "M2: t1est_deletion_with_sequence_forward_genomic",
            "M2: Specify the deleted sequence in a deletion.",
            "Switched from AL449423.14 to NG_007485.1 with other locations.",
        ],
        # "input": "AL449423.14:g.65471_65472delTC",
        "input": "NG_007485.1:g.5350_5352delCCA",
        "normalized": "NG_007485.1:g.5350_5352del",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.190_192del",
                "NG_007485.1(NP_478102.2):p.(Pro64del)",
            ),
            (
                "NG_007485.1(NM_000077.4):c.-19315_-19313del",
                "NG_007485.1(NP_000068.1):p.(=)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_with_length_forward_genomic",
            "M2: Specify the deleted sequence length in a deletion.",
            "Switched from AL449423.14 to NG_007485.1 with other locations.",
        ],
        "input": "NG_007485.1:g.5350_5352del3",
        "normalized": "NG_007485.1:g.5350_5352del",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.190_192del",
                "NG_007485.1(NP_478102.2):p.(Pro64del)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_with_sequence_reverse_coding",
            "M2: Specify the deleted sequence in a deletion on the reverse strand.",
            "Switched from AL449423.14 to NG_012337.1.",
        ],
        "input": "NG_012337.1(NM_012459.2):c.12_13delCA",
        "normalized": "NG_012337.1(NM_012459.2):c.12_13del",
        "genomic": "NG_012337.1:g.4913_4914del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_with_length_reverse_coding",
            "M2: Specify the deleted sequence length in a deletion on the reverse strand.",
            "Switched from AL449423.14 to NG_012337.1.",
        ],
        "input": "NG_012337.1(NM_012459.2):c.12_13del2",
        "normalized": "NG_012337.1(NM_012459.2):c.12_13del",
        "genomic": "NG_012337.1:g.4913_4914del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_with_sequence_reverse_ng_coding",
            "M2: Specify the deleted sequence in a deletion on the reverse strand using a genomic reference.",
        ],
        "input": "NG_008939.1:c.155_157delAAC",
        "normalized": "NG_008939.1(NM_000532.5):c.155_157del",
        "genomic": "NG_008939.1:g.5206_5208del",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.155_157del",
                "NG_008939.1(NP_000523.2):p.(Gln52del)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_with_length_reverse_ng_coding",
            "M2: Specify the deleted sequence length in a deletion on the reverse strand using a genomic reference.",
        ],
        "input": "NG_008939.1:c.155_157del3",
        "normalized": "NG_008939.1(NM_000532.5):c.155_157del",
        "genomic": "NG_008939.1:g.5206_5208del",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.155_157del",
                "NG_008939.1(NP_000523.2):p.(Gln52del)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_inversion",
            "M2: Inversion variant.",
            "Switched from AB026906.1 to NG_008939.1.",
        ],
        "input": "NG_008939.1:c.274_275inv",
        "normalized": "NG_008939.1(NM_000532.5):c.274_275inv",
        "genomic": "NG_008939.1:g.10622_10623inv",
        "coding_protein_descriptions": {
            (
                "NG_008939.1(NM_000532.5):c.274_275inv",
                "NG_008939.1(NP_000523.2):p.(Asp92Ser)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: delins_with_length",
            "Delins with explicit length of deleted sequence (bug #108).",
        ],
        "input": "NM_000193.2:c.108_109del2insG",
        "coding_protein_descriptions": {
            (
                "NM_000193.2:c.108_109delinsG",
                "NM_000193.2(NP_000184.1):p.(Lys38Serfs*2)",
            )
        },
        "to_test": True,
    },
    # test_protein_level_description: to be implemented.
    # test_protein_reference:
    # test_wnomrna_other: No equivalent reference found.
    # test_wnomrna:
    {
        "keywords": [
            "M2: test_mrna_ref_adjacent_exons_warn",
            "M2: Warning for mRNA reference where exons are not adjacent.",
            "M2: In L41870.1 exon 15 ends on 1558 and 16 starts on 1636.",
        ],
        "input": "L41870.1:c.1del",
        "normalized": "L41870.1:c.1del",
        "infos": [],  # TODO
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_mrna_ref_adjacent_exons_no_warn",
            "M2: No warning for mRNA reference where exons are adjacent.",
        ],
        "input": "NM_003002.2:c.1del",
        "normalized": "NM_003002.2:c.1del",
        "infos": [],  # TODO
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_fs_no_stop",
            "Frame shift yielding no stop codon should be described with "
            "uncertainty of the stop codon."
            "http://www.hgvs.org/mutnomen/FAQ.html#nostop",
        ],
        "input": "NM_001199.3:c.2188dup",
        "coding_protein_descriptions": {
            (
                "NM_001199.3:c.2188dup",
                "NM_001199.3(NP_001190.1):p.(Gln730Profs*?)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ext_no_stop",
            "Extension yielding no stop codon should be described with "
            "uncertainty of the stop codon."
            "http://www.hgvs.org/mutnomen/FAQ.html#nostop",
        ],
        "input": "NM_000193.2:c.1388G>C",
        "coding_protein_descriptions": {
            (
                "NM_000193.2:c.1388G>C",
                "NM_000193.2(NP_000184.1):p.(*463Serext*?)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_fs_ext_no_stop",
            "Extension yielding no stop codon should be described with "
            "uncertainty of the stop codon."
            "http://www.hgvs.org/mutnomen/FAQ.html#nostop",
            "To be implemented.",
        ],
        "input": "NM_000193.2:c.1388_1389insC",
        "coding_protein_descriptions": {
            (
                "NM_000193.2:c.1388_1389insC",
                "NM_000193.2(NP_000184.1):p.(*463Cysext*?)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_synonymous_p_is",
            "Synonymous mutation should yield a p.(=) description.",
            "Switched from AB026906.1 to NG_012337.1 and from SDHD_v001 to NM_003002.2",
        ],
        "input": "NG_012337.1(NM_003002.2):c.276C>T",
        "coding_protein_descriptions": {
            (
                # "AB026906.1(SDHD_v001):c.276C>T",
                "NG_012337.1(NM_003002.2):c.276C>T",
                "NG_012337.1(NP_002993.1):p.(=)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_synonymous_p_is_alt_start",
            "Synonymous mutation should yield a p.(=) description, also with an "
            "alternative start codon.",
        ],
        "input": "NM_024426.4:c.1107A>G",
        "coding_protein_descriptions": {
            (
                "NM_024426.4:c.1107A>G",
                "NM_024426.4(NP_077744.3):p.(=)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_start_codon",
            "Mutation of start codon should yield a p.? description.",
            "Used NG_012337.1 instead of AB026906.1." "To be implemented.",
        ],
        # "input": "AB026906.1:c.1A>G",
        "input": "NG_012337.1(NM_003002.2):c.1A>G",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_003002.2):c.1A>G",
                "NG_012337.1(NP_002993.1):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_start_codon_alt_start",
            "Mutation of start codon should yield a p.? description, also with "
            "an alternative start codon.",
            "To be implemented.",
        ],
        "input": "NM_024426.4:c.1C>G",
        "coding_protein_descriptions": {
            (
                "NM_024426.4:c.1C>G",
                "NM_024426.4(NP_077744.3):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_start_codon_yield_start_p_is",
            "Silent mutation creating new start codon should yield a p.? "
            "description. The visualisation should also render the case for "
            "the new start codon.",
            "Used NG_012337.1 instead of AB026906.1." "To be implemented.",
        ],
        # "input": "AB026906.1:c.1A>T", # yields TTG start codon
        "input": "NG_012337.1(NM_003002.2):c.1A>T",  # yields TTG start codon
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_003002.2):c.1A>T",
                "NG_012337.1(NP_002993.1):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_start_codon_alt_start_yield_start_p_is",
            "Silent mutation creating new start codon should yield a p.? "
            "description, also with an alternative start codon. The "
            "visualisation should also render the case for the new start codon.",
            "To be implemented.",
        ],
        "input": "NM_024426.4:c.1C>A",  # yields ATG start codon
        "coding_protein_descriptions": {
            (
                "NM_024426.4:c.1C>A",
                "NM_024426.4(NP_077744.3):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_start_codon_yield_start",
            "Mutation creating new start codon should yield a p.? description. "
            "The visualisation should also render the case for the new start "
            "codon.",
            "Used NG_012337.1 instead of AB026906.1." "To be implemented.",
        ],
        # "input": "AB026906.1:c.1_4delinsTTGA", # yields TTG start codon
        "input": "NG_012337.1(NM_003002.2):c.1_4delinsTTGA",
        # yields TTG start codon
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_003002.2):c.[1A>T;4G>A]",
                "NG_012337.1(NP_002993.1):p.?",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_start_codon_alt_start_yield_start",
            "Mutation creating new start codon should yield a p.? description, "
            "also with an alternative start codon. The visualisation should "
            "also render the new start codon.",
            "To be implemented.",
        ],
        "input": "NM_024426.4:c.1_4delinsATGA",  # yields ATG start codon
        "coding_protein_descriptions": {
            (
                "NM_024426.4:c.[1C>A;4C>A]",
                "NM_024426.4(NP_077744.3):p.?",
            )
        },
        "to_test": True,
    },
    # test_legend_mrna_by_construction: No longer relevant?
    {
        "keywords": [
            "M2: test_protein_ext_stop",
            "Variant in stop codon where an alternative stop codon is found "
            "downstream in the RNA should yield `ext*P` where P is a position.",
        ],
        "input": "NM_000143.3:c.1531T>G",
        "coding_protein_descriptions": {
            (
                "NM_000143.3:c.1531T>G",
                "NM_000143.3(NP_000134.2):p.(*511Glyext*3)",
            )
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_no_arg",
            "M2: Deletion without argument.",
        ],
        "input": "NM_000143.3:c.45del",
        "normalized": "NM_000143.3:c.45del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_arg_seq",
            "M2: Deletion with sequence argument.",
        ],
        "input": "NM_000143.3:c.45delG",
        "normalized": "NM_000143.3:c.45del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_arg_seq_enodna",
            "M2: Deletion with non-DNA sequence argument.",
        ],
        "input": "NM_000143.3:c.45delR",
        "to_test": False,
        "errors": ["ENODNA"],
    },
    {
        "keywords": [
            "M2: ",
            "M2: ",
        ],
        "input": "",
        "normalized": "",
        "to_test": False,
    },
]
