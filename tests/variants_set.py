M2_TESTS = [
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
                "NG_007485.1(NM_058195.3):c.204_206del",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGlu)",
            ),
        },
        "protein_description": "NG_007485.1(NP_000068.1):p.(Met54_Gly55delinsSer)",
        # "AL449423.14(CDKN2A_i001):p.(Met54_Gly55delinsSer)",
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
                "NG_007485.1(NM_058195.3):c.204_205insATC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69insIle)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444insGAT"],
        "protein_description": "NG_007485.1(NP_000068.1):p.(Met54delinsIleSer)",
        # "AL449423.14(CDKN2A_i001):p.(Met54delinsIleSer)"
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
                "NG_007485.1(NM_058195.3):c.204_205insATC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69insIle)",
            ),
        },
        "protein_description": "NG_007485.1(NP_000068.1):p.(Met54delinsIleSer)",
        # "AL449423.14(CDKN2A_i001):p.(Met54delinsIleSer)"
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
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444delinsGGGAT"],
        "protein_description": "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
        # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
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
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444delinsGGGAT"],
        "protein_description": "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
        # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
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
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444delinsGGGAT"],
        "protein_description": "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
        # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
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
                "NG_007485.1(NM_058195.3):c.204_205delinsATCCC",
                "NG_007485.1(NP_478102.2):p.(Asp68_Gly69delinsGluSerArg)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+3443_616+3444delinsGGGAT"],
        "protein_description": "NG_007485.1(NP_000068.1):p.(Met54delinsAsnPro)",
        # "AL449423.14(CDKN2A_i001):p.(Met54delinsAsnPro)"
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
        "infos": [
            "ICORRECTEDCOORDINATESYSTEM",
            "IWHOLETRANSCRIPTEXON",
            "IMRNAGENOMICTIP",
        ],
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
        "protein_description": "NM_003002.2(NP_002993.1):p.(Asp92Thrfs*43)",
        "infos": ["IWHOLETRANSCRIPTEXON", "IMRNAGENOMICTIP"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_no_roll",
            "M2: Just a variant where we cannot roll.",
        ],
        "input": "NM_003002.2:c.274del",
        "normalized": "NM_003002.2:c.274del",
        "protein_description": "NM_003002.2(NP_002993.1):p.(Asp92Thrfs*43)",
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
        "normalized": "NM_000088.3:n.460del",
        "to_test": False,
    },
    {
        "keywords": [
            "M2: test_roll_both_ins",
            """M2:
            Insertion that rolls should not use the same inserted sequence in
            descriptions on forward and reverse strands.
            Here we have the following situation on the forward strand:
                                      4812 (genomic)
                                        |
                CCGCACCTGTGCAGTAAACTGCGCCTTCTGCTGCTCGGCGGCCACCAGGC
            Now, an insertion of TAC after 65470 should be rolled to an insertion
            of ACT after 65471:

                CCGCACCTGTGCAGTAAACTGCGCC --- TTCTGCTGCTCGGCGGCCACCAGGC
                CCGCACCTGTGCAGTAAACTGCGCC TAC TTCTGCTGCTCGGCGGCCACCAGGC  =>

                CCCGCACCTGTGCAGTAAACTGCGC --- CTTCTGCTGCTCGGCGGCCACCAGG
                CCCGCACCTGTGCAGTAAACTGCGC CTA CTTCTGCTGCTCGGCGGCCACCAGG
            However, in NM_012459.2 (on the reverse strand), this insertion should
            roll the other direction and the inserted sequence should be the reverse
            complement of CTA, which is TAG, and not that of ACT, which is AGT.
            The next test (test_roll_reverse_ins) tests the situation for an input
            of AL449423.14:g.65471_65472insACT, where only the reverse roll should
            be done.
            """,
            "Switched from AL449423.14 to NG_012337.1 and updated locations.",
        ],
        # "input": "AL449423.14:g.65470_65471insTAC",
        "input": "NG_012337.1:g.4812_4813insTAC",
        "normalized": "NG_012337.1:g.4813_4814insACT",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_018195.3):c.*3687_*3688insACT",
                "NG_012337.1(NP_060665.3):p.(=)",
            ),
            (
                "NG_012337.1(NM_001082969.1):c.*3687_*3688insACT",
                "NG_012337.1(NP_001076438.1):p.(=)",
            ),
            (
                "NG_012337.1(NM_001082970.1):c.*3687_*3688insACT",
                "NG_012337.1(NP_001076439.1):p.(=)",
            ),
            (
                "NG_012337.1(NM_001082970.1):c.*3687_*3688insACT",
                "NG_012337.1(NP_001076439.1):p.(=)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_roll_reverse_ins",
            """M2:
            Insertion that rolls on the reverse strand should not use the same
            inserted sequence in descriptions on forward and reverse strands.
            """,
            "Switched from AL449423.14 to NG_007485.1 and updated locations."
            "Note that the noncoding description is on the reverse strand",
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
        "to_test": True,
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
        "protein_description": "NM_000143.3(NP_000134.2):p.(=)",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_ins_cds_start_after",
            "M2: Insertion after CDS start boundary should be included in CDS.",
        ],
        "input": "NM_000143.3:c.1_2insCAT",
        "normalized": "NM_000143.3:c.1_2insCAT",
        "protein_description": "NM_000143.3(NP_000134.2):p.?",
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
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_exon",
            "M2: Deletion of an entire exon should be possible.",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.632-5_681+7del",
        "normalized": "NG_012772.1(NM_000059.3):c.632-5_681+7del",
        "genomic": "NG_012772.1:g.18959_19020del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_exon_exact",
            "M2: Deletion of exactly an exon should be possible.",
        ],
        "input": "NG_012772.1(BRCA2_v001):c.632_681del",
        "normalized": "NG_012772.1(NM_000059.3):c.632_681del",
        "protein_description": "NG_012772.1(NP_000050.2):p.(Val211Glufs*10)",
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
        "protein_description": "NG_012772.1(NP_000050.2):p.(Val208Tyrfs*3)",
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
        "protein_description": "NG_012772.1(NP_000050.2):p.(Val208_Asp224del)",
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
        "protein_description": "NG_012772.1(NP_000050.2):p.(Val211Glufs*10)",
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
        "protein_description": "NG_012772.1(NP_000050.2):p.(Asp23_Leu105del)",
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
        "protein_description": "NM_000143.3(NP_000134.2):p.(Glu247Alafs*27)",
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
        "protein_description": "NG_012337.1(NP_036591.2):p.(Ser5Aspfs*21)",
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
        "protein_description": "NG_012337.1(NP_036591.2):p.(Ser5Aspfs*21)",
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
        "protein_description": "NG_012337.1(NP_036591.2):p.(Ser5Phefs*22)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Gln52_Arg53insValLeuCysSerLeuSerGly)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Arg53Alafs*5)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Arg53Phefs*11)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Arg53Alafs*5)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Arg53Phefs*11)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Arg53_Arg54delinsSerCysAlaHisTyrLeuAla)",
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
        "protein_description": "LRG_1(p1):p.(Gly89Val)",
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
            "Note: NM_000532.5 is not on the reverse strand.",
        ],
        "input": "NG_008939.1:c.155_157delAAC",
        "normalized": "NG_008939.1(NM_000532.5):c.155_157del",
        "genomic": "NG_008939.1:g.5206_5208del",
        "protein_description": "NG_008939.1(NP_000523.2):p.(Gln52del)",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_deletion_with_length_reverse_ng_coding",
            "M2: Specify the deleted sequence length in a deletion on the reverse strand using a genomic reference.",
            "Note: NM_000532.5 is not on the reverse strand.",
        ],
        "input": "NG_008939.1:c.155_157del3",
        "normalized": "NG_008939.1(NM_000532.5):c.155_157del",
        "genomic": "NG_008939.1:g.5206_5208del",
        "protein_description": "NG_008939.1(NP_000523.2):p.(Gln52del)",
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
        "protein_description": "NG_008939.1(NP_000523.2):p.(Asp92Ser)",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: delins_with_length",
            "Delins with explicit length of deleted sequence (bug #108).",
        ],
        "input": "NM_000193.2:c.108_109del2insG",
        "protein_description": "NM_000193.2(NP_000184.1):p.(Lys38Serfs*2)",
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
        "protein_description": "NM_001199.3(NP_001190.1):p.(Gln730Profs*?)",
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
        "protein_description": "NM_000193.2(NP_000184.1):p.(*463Serext*?)",
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
        "protein_description": "NM_000193.2(NP_000184.1):p.(*463Cysext*?)",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_synonymous_p_is",
            "Synonymous mutation should yield a p.(=) description.",
            "Switched from AB026906.1 to NG_012337.1 and from SDHD_v001 to NM_003002.2",
        ],
        "input": "NG_012337.1(NM_003002.2):c.276C>T",
        # "AB026906.1(SDHD_v001):c.276C>T"
        "protein_description": "NG_012337.1(NP_002993.1):p.(=)",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_synonymous_p_is_alt_start",
            "Synonymous mutation should yield a p.(=) description, also with an "
            "alternative start codon.",
        ],
        "input": "NM_024426.4:c.1107A>G",
        "protein_description": "NM_024426.4(NP_077744.3):p.(=)",
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
        "protein_description": "NG_012337.1(NP_002993.1):p.?",
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
        "protein_description": "NM_024426.4(NP_077744.3):p.?",
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
        "protein_description": "NG_012337.1(NP_002993.1):p.?",
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
        "protein_description": "NM_024426.4(NP_077744.3):p.?",
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
        "protein_description": "NG_012337.1(NP_002993.1):p.?",
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
        "protein_description": "NM_024426.4(NP_077744.3):p.?",
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
        "protein_description": "NM_000143.3(NP_000134.2):p.(*511Glyext*3)",
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
        "input": "NM_000143.3:c.45delU",
        "errors": ["ENODNA"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_arg_seq_eref",
            "M2: Deletion with non-reference sequence argument.",
        ],
        "input": "NM_000143.3:c.45delT",
        "errors": ["ESEQUENCEMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_arg_len",
            "M2: Deletion with length argument.",
        ],
        "input": "NM_000143.3:c.45del1",
        "normalized": "NM_000143.3:c.45del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_arg_len_earglen",
            "M2: Deletion with incorrect length argument.",
        ],
        "input": "NM_000143.3:c.45del4",
        "errors": ["ELENGTHMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_range_no_arg",
            "M2: Range deletion without argument.",
        ],
        "input": "NM_000143.3:c.44_47del",
        "normalized": "NM_000143.3:c.44_47del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_range_arg_seq",
            "M2: Range deletion with sequence argument.",
        ],
        "input": "NM_000143.3:c.44_47delTGCG",
        "normalized": "NM_000143.3:c.44_47del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_range_arg_seq_enodna",
            "M2: Range deletion with non-DNA sequence argument.",
        ],
        "input": "NM_000143.3:c.44_47delTCUG",
        "errors": ["ENODNA"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_range_arg_seq_eref",
            "M2: Range deletion with non-reference sequence argument.",
        ],
        "input": "NM_000143.3:c.44_47delTGGG",
        "errors": ["ESEQUENCEMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_range_arg_len",
            "M2: Range deletion with length argument.",
        ],
        "input": "NM_000143.3:c.44_47del4",
        "normalized": "NM_000143.3:c.44_47del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_range_arg_len_earglen",
            "M2: Range deletion with incorrect length argument.",
        ],
        "input": "NM_000143.3:c.44_47del3",
        "errors": ["ELENGTHMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_no_arg",
            "M2: Delins without argument.",
        ],
        "input": "NM_000143.3:c.45delinsATC",
        "normalized": "NM_000143.3:c.45delinsATC",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_arg_seq",
            "M2: Delins with sequence argument.",
        ],
        "input": "NM_000143.3:c.45delGinsATC",
        "normalized": "NM_000143.3:c.45delinsATC",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_arg_seq_enodna",
            "M2: Delins with non-DNA sequence argument.",
        ],
        "input": "NM_000143.3:c.45deluinsATC",
        "errors": ["ENODNA"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_arg_seq_eref",
            "M2: Delins with non-reference sequence argument.",
        ],
        "input": "NM_000143.3:c.45delTinsATC",
        "errors": ["ESEQUENCEMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_arg_len",
            "M2: Delins with length argument.",
        ],
        "input": "NM_000143.3:c.45del1insATC",
        "normalized": "NM_000143.3:c.45delinsATC",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_arg_len_earglen",
            "M2: Delins with incorrect length argument.",
        ],
        "input": "NM_000143.3:c.45del4insATC",
        "errors": ["ELENGTHMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_no_arg",
            "M2: Range delins without argument.",
        ],
        "input": "NM_000143.3:c.44_47delinsATC",
        "normalized": "NM_000143.3:c.44_47delinsATC",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_arg_seq",
            "M2: Range delins with sequence argument.",
        ],
        "input": "NM_000143.3:c.44_47delTGCGinsATC",
        "normalized": "NM_000143.3:c.44_47delinsATC",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_arg_seq_enodna",
            "M2: Range delins with non-DNA sequence argument.",
        ],
        "input": "NM_000143.3:c.44_47delTCuGinsATC",
        "errors": ["ENODNA"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_arg_seq_eref",
            "M2: Range delins with non-reference sequence argument.",
        ],
        "input": "NM_000143.3:c.44_47delTGGGinsATC",
        "errors": ["ESEQUENCEMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_arg_len",
            "M2: Range delins with length argument.",
        ],
        "input": "NM_000143.3:c.44_47del4insATC",
        "normalized": "NM_000143.3:c.44_47delinsATC",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_arg_len_earglen",
            "M2: Range delins with incorrect length argument.",
        ],
        "input": "NM_000143.3:c.44_47del3insATC",
        "errors": ["ELENGTHMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_no_arg",
            "M2: Duplication without argument.",
        ],
        "input": "NM_000143.3:c.45dup",
        "normalized": "NM_000143.3:c.45dup",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_arg_seq",
            "M2: Duplication with sequence argument.",
        ],
        "input": "NM_000143.3:c.45dupG",
        "normalized": "NM_000143.3:c.45dup",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_arg_seq_enodna",
            "M2: Duplication with non-DNA sequence argument.",
        ],
        "input": "NM_000143.3:c.45dupu",
        "errors": ["ENODNA"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_arg_seq_eref",
            "M2: Duplication with non-reference sequence argument.",
        ],
        "input": "NM_000143.3:c.45dupT",
        "errors": ["ESEQUENCEMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_arg_len",
            "M2: Duplication with length argument.",
        ],
        "input": "NM_000143.3:c.45dup1",
        "normalized": "NM_000143.3:c.45dup",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_arg_len_earglen",
            "M2: Duplication with incorrect length argument.",
        ],
        "input": "NM_000143.3:c.45dup4",
        "errors": ["ELENGTHMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_range_no_arg",
            "M2: Range duplication without argument.",
        ],
        "input": "NM_000143.3:c.44_47dup",
        "normalized": "NM_000143.3:c.44_47dup",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_range_arg_seq",
            "M2: Range duplication with sequence argument.",
        ],
        "input": "NM_000143.3:c.44_47dupTGCG",
        "normalized": "NM_000143.3:c.44_47dup",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_range_arg_seq_enodna",
            "M2: Range duplication with non-DNA sequence argument.",
        ],
        "input": "NM_000143.3:c.44_47dupTCuG",
        "errors": ["ENODNA"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_range_arg_seq_eref",
            "M2: Range duplication with non-reference sequence argument.",
        ],
        "input": "NM_000143.3:c.44_47dupTGGG",
        "errors": ["ESEQUENCEMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_range_arg_len",
            "M2: Range duplication with length argument.",
        ],
        "input": "NM_000143.3:c.44_47dup4",
        "normalized": "NM_000143.3:c.44_47dup",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_range_arg_len_earglen",
            "M2: Range duplication with incorrect length argument.",
        ],
        "input": "NM_000143.3:c.44_47dup3",
        "errors": ["ELENGTHMISMATCH"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_arg_seq_reverse",
            "M2: Deletion with sequence argument (reverse strand).",
            "Switched from AL449423.14 to NG_012337.1",
        ],
        "input": "NG_012337.1(TIMM8B):c.12delC",
        "normalized": "NG_012337.1(NM_012459.2):c.12del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_del_range_arg_seq_reverse(",
            "M2: Range deletion with sequence argument (reverse strand).",
        ],
        "input": "NG_012337.1(TIMM8B):c.12_15delCAGC",
        "normalized": "NG_012337.1(NM_012459.2):c.12_15del",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_arg_seq_reverse",
            "M2: Delins with sequence argument (reverse strand).",
        ],
        "input": "NG_012337.1(TIMM8B):c.12delCinsAT",
        "normalized": "NG_012337.1(NM_012459.2):c.12delinsAT",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_delins_range_arg_seq_reverse",
            "M2: Range delins with sequence argument (reverse strand).",
        ],
        "input": "NG_012337.1(TIMM8B):c.12_15delCAGCinsTTT",
        "normalized": "NG_012337.1(NM_012459.2):c.12_15delinsTTT",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_arg_seq_reverse",
            "M2: Duplication with sequence argument (reverse strand).",
        ],
        "input": "NG_012337.1(TIMM8B):c.12dupC",
        "normalized": "NG_012337.1(NM_012459.2):c.12dup",
        "genomic": "NG_012337.1:g.4911dup",
        "coding_protein_descriptions": {
            ("NG_012337.1(NM_018195.3):c.*3785dup", "NG_012337.1(NP_060665.3):p.(=)"),
            (
                "NG_012337.1(NM_001082969.1):c.*3785dup",
                "NG_012337.1(NP_001076438.1):p.(=)",
            ),
            (
                "NG_012337.1(NM_001082970.1):c.*3785dup",
                "NG_012337.1(NP_001076439.1):p.(=)",
            ),
        },
        "protein_description": "NG_012337.1(NP_036591.2):p.(His4Glnfs*21)",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_dup_range_arg_seq_reverse",
            "M2: Range duplication with sequence argument (reverse strand).",
        ],
        "input": "NG_012337.1(TIMM8B):c.12_15dupCAGC",
        "normalized": "NG_012337.1(NM_012459.2):c.12_15dup",
        "genomic": "NG_012337.1:g.4908_4911dup",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_018195.3):c.*3782_*3785dup",
                "NG_012337.1(NP_060665.3):p.(=)",
            ),
            (
                "NG_012337.1(NM_001082969.1):c.*3782_*3785dup",
                "NG_012337.1(NP_001076438.1):p.(=)",
            ),
            (
                "NG_012337.1(NM_001082970.1):c.*3782_*3785dup",
                "NG_012337.1(NP_001076439.1):p.(=)",
            ),
        },
        "protein_description": "NG_012337.1(NP_036591.2):p.(His4Glnfs*22)",
        "to_test": True,
    },
    {
        "keywords": [
            "M2: test_accno_as_transcript_variant",
            "M2: Test accession number as transcript variant identifier is accepted.",
        ],
        "input": "NG_012337.1(NM_012459.2):c.12_13insGATC",
        "normalized": "NG_012337.1(NM_012459.2):c.12_13insGATC",
        "genomic": "NG_012337.1:g.4911_4912insATCG",
        "to_test": True,
    },
    {
        "keywords": ["M2: `test_accno_as_transcript_variant` the error part"],
        "input": "NG_012337.1(DUMMYACCNO_9999.9):c.12_13insGATC",
        "errors": ["ENOSELECTORFOUND"],
        "to_test": True,
    },
    {
        "keywords": ["protein", "M2: similar to `accno_as_transcript_variant`"],
        "input": "NG_012337.1(NM_012459.2):c.4_5insGTA",
        "protein_description": "NG_012337.1(NP_036591.2):p.(Arg2_Lys3insSer)",
        "to_test": True,
    },
    {
        "keywords": ["protein"],
        "input": "LRG_199t1:c.235_237delinsTAT",
        "genomic": "LRG_199:g.[499798A>T;499800G>T]",
        "normalized": "LRG_199(t1):c.[235A>T;237G>T]",
        "protein_description": "LRG_199(p1):p.(Lys79Tyr)",
        "to_test": True,
    },
    # All the IVS/EX tests are not any longer relevant.
]


TESTS = [
    {
        "keywords": [
            "correct input",
            "no errors",
            "no warnings",
            "genomic reference",
            "NCBI",
            "g.",
            "single variant",
            "substitution",
        ],
        "input": "NG_012337.1:g.4C>T",
        "normalized": "NG_012337.1:g.4C>T",
        "to_test": True,
    },
    {
        "keywords": [
            "correct input",
            "no errors",
            "no warnings",
            "genomic reference",
            "NCBI",
            "g.",
            "single variant",
            "deletion range",
        ],
        "input": "NG_017013.2:g.17013_17014del",
        "normalized": "NG_017013.2:g.17013_17014del",
        "to_test": True,
    },
    {
        "keywords": [
            "correct input",
            "no errors",
            "no warnings",
            "genomic reference",
            "NCBI",
            "g.",
            "single variant",
            "deletion insertion single range",
        ],
        "input": "NG_012337.1:g.4delins7_31",
        "normalized": "NG_012337.1:g.4delins7_31",
        "to_test": True,
    },
    {
        "keywords": [
            "correct input",
            "no errors",
            "no warnings",
            "genomic reference",
            "NCBI",
            "g.",
            "single variant",
            "duplication single",
        ],
        "input": "NG_017013.2:g.19258dup",
        "normalized": "NG_017013.2:g.19258dup",
        "to_test": True,
    },
    {
        "keywords": [
            "correct input",
            "no errors",
            "no warnings",
            "genomic reference",
            "NCBI",
            "g.",
            "single variant",
            "duplication range",
        ],
        "input": "NG_017013.2:g.16508_16509dup",
        "normalized": "NG_017013.2:g.16508_16509dup",
        "to_test": True,
    },
    {
        "keywords": [
            "correct input",
            "no errors",
            "no warnings",
            "genomic reference",
            "NCBI",
            "g.",
            "single variant",
            "insertion single sequence",
        ],
        "input": "NG_012337.1:g.90_91insC",
        "normalized": "NG_012337.1:g.90_91insC",
        "to_test": True,
    },
    # ---------
    {
        "keywords": [
            "correct input",
            "no errors",
            "no warnings",
            "NCBI",
            "c.",
            "single variant",
            "substitution",
            "roll forward",
        ],
        "input": "NG_007485.1(NM_000077.4):c.274G>T",
        "normalized": "NG_007485.1(NM_000077.4):c.274G>T",
        "genomic": "NG_007485.1:g.28407G>T",
        "to_test": True,
    },
    # ---------
    {
        "keywords": [
            "corrected",
            "no errors",
            "warnings",
            "NCBI",
            "g.",
            "single variant",
            "deletion range",
            "roll forward",
        ],
        "input": "NG_012337.1:g.26_31del",
        "normalized": "NG_012337.1:g.29_34del",
        "warnings": "Deletion at position 26_31 was given, however, the HGVS "
        "notation prescribes that on the forward strand it should "
        "be at position 29_34.",
        "to_test": True,
    },
    {
        "keywords": [
            "corrected",
            "no errors",
            "warnings",
            "NCBI",
            "g.",
            "single variant",
            "deletion range",
            "roll forward",
        ],
        "input": "NG_017013.2:g.17011_17012del",
        "normalized": "NG_017013.2:g.17013_17014del",
        "to_test": True,
    },
    # ---------
    {
        "keywords": [],
        "input": "NG_017013.2:g.17415_17417delinsGCG",
        "normalized": "NG_017013.2:g.[17415C>G;17417A>G]",
        "genomic": "NG_017013.2:g.[17415C>G;17417A>G]",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.17496_17497insAGCTGCTCAGATAGCGA",
        "normalized": "NG_017013.2:g.17496_17497ins[A;17481_17496]",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.1_32772del",
        "normalized": "NG_017013.2:g.1_32772del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.[16985A>T;17013_17014del]",
        "normalized": "NG_017013.2:g.[16985A>T;17013_17014del]",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.[16985A>T;17011_17012del]",
        "normalized": "NG_017013.2:g.[16985A>T;17013_17014del]",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.17011_17012del",
        "normalized": "NG_017013.2:g.17013_17014del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.4delins7_50",
        "normalized": "NG_012337.1:g.[3_4insGGTT;5_6insACCATATCTCTACTTTGTGTTTATGTTTGTGTATGCATT]",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.26_31del",
        "normalized": "NG_012337.1:g.29_34del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.100_200delins100_101",
        "normalized": "NG_012337.1:g.102_200del",
        "genomic": "NG_012337.1:g.102_200del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.17471_17471del",
        "normalized": "NG_017013.2:g.17471del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.18748_18750delinsCAT",
        "normalized": "NG_017013.2:g.18749G>A",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.17394_17395insC",
        "normalized": "NG_017013.2:g.17394dup",
        "genomic": "NG_017013.2:g.17394dup",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.274G>T",
        "normalized": "NG_012337.1(NM_003002.2):c.274G>T",
        "genomic": "NG_012337.1:g.7125G>T",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.1del",
        "normalized": "NG_012337.1(NM_003002.2):c.1del",
        "genomic": "NG_012337.1:g.5062del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.-1del",
        "normalized": "NG_012337.1(NM_003002.2):c.-1del",
        "genomic": "NG_012337.1:g.5061del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.52+1del",
        "normalized": "NG_012337.1(NM_003002.2):c.52+1del",
        "genomic": "NG_012337.1:g.5114del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.*824del",
        "normalized": "NG_012337.1(NM_003002.2):c.*824del",
        "genomic": "NG_012337.1:g.13948del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.*824+10del",
        "normalized": "NG_012337.1(NM_003002.2):c.*834del",
        "genomic": "NG_012337.1:g.13958del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NM_003002.4:c.1del",
        "normalized": "NM_003002.4:c.1del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_029724.1(NM_004321.7):c.101del",
        "normalized": "NG_029724.1(NM_004321.7):c.102del",
        "genomic": "NG_029724.1:g.27557del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_029724.1(NM_004321.7):c.101del",
        "normalized": "",
        "genomic": "NG_029724.1:g.27557del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_029724.1:g.10_20del",
        "normalized": "NG_029724.1:g.10_20del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_029724.1:g.10_20del11",
        "normalized": "NG_029724.1:g.10_20del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_029724.1:g.10delG",
        "normalized": "NG_029724.1:g.10del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_029724.1:g.10_11delGT",
        "normalized": "NG_029724.1:g.10_11del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1(NM_022124.6):c.3306del",
        "normalized": "NG_008835.1(NM_022124.6):c.3306del",
        "genomic": "NG_008835.1:g.320804del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1(NM_022124.6):c.3304del",
        "normalized": "NG_008835.1(NM_022124.6):c.3306del",
        "genomic": "NG_008835.1:g.320804del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1(NM_022124.6):c.3304_3305del",
        "normalized": "NG_008835.1(NM_022124.6):c.3305_3306del",
        "genomic": "NG_008835.1:g.320803_320804del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1(NM_001168390.2):c.*3188del",
        "normalized": "NG_008835.1(NM_001168390.2):c.*3188del",
        "genomic": "NG_008835.1:g.320804del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1(NM_001168390.2):c.*3187del",
        "normalized": "NG_008835.1(NM_001168390.2):c.*3188del",
        "genomic": "NG_008835.1:g.320804del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1(NM_001168390.2):c.*3186del",
        "normalized": "NG_008835.1(NM_001168390.2):c.*3188del",
        "genomic": "NG_008835.1:g.320804del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1(NM_001168390.2):c.*3186_*3187del",
        "normalized": "NG_008835.1(NM_001168390.2):c.*3187_*3188del",
        "genomic": "NG_008835.1:g.320803_320804del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1:g.320801del",
        "normalized": "NG_008835.1:g.320801del",
        "genomic": "NG_008835.1:g.320801del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1:g.320802del",
        "normalized": "NG_008835.1:g.320804del",
        "genomic": "NG_008835.1:g.320804del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1:g.320801_320802del",
        "normalized": "NG_008835.1:g.320801_320802del",
        "genomic": "NG_008835.1:g.320801_320802del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1:g.320802_320803del",
        "normalized": "NG_008835.1:g.320803_320804del",
        "genomic": "NG_008835.1:g.320803_320804del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_008835.1(CDH23_v001):c.1449+846delA",
        "normalized": "NG_008835.1(NM_022124.6):c.1449+858del",
        "genomic": "NG_008835.1:g.255529del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_009113.2(NM_014249.4):c.948delC",
        "normalized": "NG_009113.2(NM_014249.4):c.951del",
        "genomic": "NG_009113.2:g.8038del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_007107.2(NM_004992.3):c.378-17delT",
        "normalized": "NG_007107.2(NM_004992.3):c.378-17del",
        "genomic": "NG_007107.2:g.110661del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_009113.2(NR2E3_v002):c.948delC",
        "normalized": "",
        "genomic": "NG_009113.2:g.8038del",
        "to_test": False,
    },
    {
        "keywords": [],
        "input": "NG_009113.2(NM_016346.4):c.948delC",
        "normalized": "NG_009113.2(NM_016346.4):c.951del",
        "genomic": "NG_009113.2:g.8038del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_009497.1(NM_206933.2):c.8682-19dupT",
        "normalized": "NG_009497.1(NM_206933.2):c.8682-19dup",
        "genomic": "NG_009497.1:g.561208dup",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_009497.1(USH2A_v001):c.8682-19dup",
        "normalized": "NG_009497.1(NM_206933.2):c.8682-19dup",
        "genomic": "NG_009497.1:g.561208dup",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_007485.1(NM_058195.3):c.141_142del",
        "normalized": "NG_007485.1(NM_058195.3):c.141_142del",
        "genomic": "NG_007485.1:g.5301_5302del",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_000077.4):c.-19364_-19363del",
                "NG_007485.1(NP_000068.1):p.(=)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+26436_616+26437del"],
        "protein_description": "NG_007485.1(NP_478102.2):p.(Met48Alafs*14)",
        "rna_description": "NG_007485.1(NM_058195.3):r.(141_142del)",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.169T>A",
        "normalized": "NG_012337.1(NM_003002.2):c.169T>A",
        "genomic": "NG_012337.1:g.6127T>A",
        "protein_description": "NG_012337.1(NP_002993.1):p.(Ser57Thr)",
        "rna_description": "NG_012337.1(NM_003002.2):r.(169u>a)",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6932_6933insAGCAACGTGATCGCCTCCCTCACCT",
        "normalized": "LRG_303:g.6908_6932dup",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6932_6933insAGCAACGTGATCGCCTCCCTCACCT[0]",
        "normalized": "LRG_303:g.=",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6932_6933insAGCAACGTGATCGCCTCCCTCACCT[1]",
        "normalized": "LRG_303:g.6908_6932dup",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6932_6933insAGCAACGTGATCGCCTCCCTCACCT[2]",
        "normalized": "LRG_303:g.6908_6932AGCAACGTGATCGCCTCCCTCACCT[3]",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6908_6932AGCAACGTGATCGCCTCCCTCACCT[0]",
        "normalized": "LRG_303:g.6908_6932del",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.4_5insGT[5]",
        "normalized": "LRG_303:g.1_4GT[7]",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.3_4GT[5]",
        "normalized": "LRG_303:g.1_4GT[6]",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.1_2GT[0]",
        "normalized": "LRG_303:g.3_4del",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.3_4GT[0]",
        "normalized": "LRG_303:g.3_4del",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.3_4GT[2]",
        "normalized": "LRG_303:g.3_4dup",
        "normalized-alt": "LRG_303:g.1_4[3]",  # TODO: check the alt.
        "to_test": True,
    },
    {
        "keywords": ["repeat", "reverse_strand"],
        "input": "NG_009299.1(NM_017668.3):c.*16_*18T[8]",
        "normalized": "NG_009299.1(NM_017668.3):c.*16_*19T[9]",
        "genomic": "NG_009299.1:g.137761_137764A[9]",
        "to_test": True,
    },
    {
        "keywords": ["repeat", "reverse_strand"],
        "input": "NG_009299.1(NM_017668.3):c.33_35CAA[5]",
        "normalized": "NG_009299.1(NM_017668.3):c.34_39AAC[6]",
        "genomic": "NG_009299.1:g.137803_137808TTG[6]",
        "to_test": True,
    },
    {
        "keywords": ["issue #25"],
        "input": "NG_012337.1:g.[7109T>A;7110del]",
        "normalized": "NG_012337.1:g.7109_7110delinsA",
        "genomic": "NG_012337.1:g.7109_7110delinsA",
        "to_test": True,
    },
    {
        "keywords": ["issue #10"],
        "input": "NG_012337.1(SDHD_v001):c.274G>T",
        "normalized": "NG_012337.1(NM_003002.2):c.274G>T",
        "genomic": "NG_012337.1:g.7125G>T",
        "to_test": True,
    },
    {
        "keywords": ["issue #10"],
        "input": "NG_012337.1(SDHD):c.274G>T",
        "normalized": "NG_012337.1(NM_003002.2):c.274G>T",
        "genomic": "NG_012337.1:g.7125G>T",
        "to_test": True,
    },
    {
        "keywords": ["issue #10"],
        "input": "NG_012337.1(10683):c.274G>T",
        "normalized": "NG_012337.1(NM_003002.2):c.274G>T",
        "genomic": "NG_012337.1:g.7125G>T",
        "to_test": True,
    },
    {
        "keywords": ["reverse strand sequence"],
        "input": "NG_009299.1(NM_017668.3):c.41A>C",
        "normalized": "NG_009299.1(NM_017668.3):c.41A>C",
        "genomic": "NG_009299.1:g.137800T>G",
        "coding_protein_descriptions": {
            (
                "NG_009299.1(NM_001040113.2):c.4316T>G",
                "NG_009299.1(NP_001035202.1):p.(Leu1439Arg)",
            ),
            (
                "NG_009299.1(NM_002474.3):c.4295T>G",
                "NG_009299.1(NP_002465.1):p.(Leu1432Arg)",
            ),
            (
                "NG_009299.1(NM_001143979.2):c.41A>C",
                "NG_009299.1(NP_001137451.1):p.(Gln14Profs*68)",
            ),
        },
        "protein_description": "NG_009299.1(NP_060138.1):p.(Gln14Profs*68)",
        "rna_description": "NG_009299.1(NM_017668.3):r.(40a>c)",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.7125G>T",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_003002.2):c.274G>T",
                "NG_012337.1(NP_002993.1):p.(Asp92Tyr)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "LRG_24:g.5526_5533del",
        "coding_protein_descriptions": {
            (
                "LRG_24(t1):c.127_134del",
                "LRG_24(p1):p.(Gly43Argfs*65)",
            ),
            (
                "LRG_24(t2):c.127_134del",
                "LRG_24(p2):p.(Gly43Argfs*65)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_012459.2):c.5_6delinsTAG",
        "protein_description": "NG_012337.1(NP_036591.2):p.(Arg2Leufs*23)",
        "rna_description": "NG_012337.1(NM_012459.2):r.(5_6delinsuag)",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_012459.2):c.4_6delinsGTA",
        "protein_description": "NG_012337.1(NP_036591.2):p.(Arg2Val)",
        "rna_description": "NG_012337.1(NM_012459.2):r.(4_6delinsgua)",
        "to_test": True,
    },
    {
        "keywords": ["reference", "LRG", "replace"],
        "input": "LRG_303(t1):c.10_11insLRG_1t1:c.100_101",
        "normalized": "LRG_303(t1):c.10_11insGA",
        "infos": ["ICORRECTEDLRGREFERENCE"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "LRG_303:g.[105_106del;6681G>C;6883_6884insTTTCGCCCCTTTCGCCCC]",
        "normalized": "LRG_303:g.[108_109del;6681G>C;6875_6883TTTCGCCCC[3]]",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.13124_13125del",
        "normalized": "NG_012337.1:g.13124_13125del",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_003002.2):c.480_*1del",
                "NG_012337.1(NP_002993.1):p.(*160Cysext*29)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": ["reverse strand"],
        "input": "NG_012337.1(NM_012459.2):c.[10del;23_25del;36_37insAAT]",
        "normalized": "NG_012337.1(NM_012459.2):c.[10del;23_25del;37_38insATA]",
        "genomic": "NG_012337.1:g.[4886_4887insATT;4898_4900del;4913del]",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.7125G>TA",
        "normalized": "NG_012337.1:g.7125delinsTA",
        "infos": ["ICORRECTEDVARIANTTYPE"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.[274G>T;278A>G]",
        "normalized": "NG_012337.1(NM_003002.2):c.[274G>T;278A>G]",
        "genomic": "NG_012337.1:g.[7125G>T;7129A>G]",
        "protein_description": "NG_012337.1(NP_002993.1):p.(Asp92_Tyr93delinsTyrCys)",
        "rna_description": "NG_012337.1(NM_003002.2):r.([274g>u;278a>g])",
        "to_test": True,
    },
    {
        "keywords": ["no operation"],
        "input": "NG_012337.1:274",
        "normalized": "NG_012337.1:g.274",
        "to_test": True,
    },
    {
        "keywords": ["no operation"],
        "input": "NG_012337.1(NM_003002.2):c.[274;600]",
        "normalized": "NG_012337.1(NM_003002.2):c.[274;*120]",
        "genomic": "NG_012337.1:g.[7125;13244]",
        "to_test": True,
    },
    {
        "keywords": ["protein reverse strand end of CDS"],
        "input": "NG_012337.1(NM_012459.2):c.297_*1del",
        "normalized": "NG_012337.1(NM_012459.2):c.297_*1del",
        "genomic": "NG_012337.1:g.3448_3449del",
        "protein_description": "NG_012337.1(NP_036591.2):p.(*99Tyrext*6)",
        "rna_description": "NG_012337.1(NM_012459.2):r.(297_*1del)",
        "to_test": True,
    },
    {
        "keywords": ["protein reverse strand end of CDS"],
        "input": "NG_012337.1(NM_012459.2):c.-35_*1del",
        "normalized": "NG_012337.1(NM_012459.2):c.-35_*1del",
        "genomic": "NG_012337.1:g.3448_4957del",
        "to_test": True,
    },
    {
        "keywords": ["protein reverse strand end of CDS"],
        "input": "NG_012337.1(NM_012459.2):c.-1_*1del",
        "normalized": "NG_012337.1(NM_012459.2):c.1_*2del",
        "genomic": "NG_012337.1:g.3449_4924del",
        "protein_description": "NG_012337.1(NP_036591.2):p.?",
        "to_test": True,
    },
    {
        "keywords": ["reverse strand"],
        "input": "NG_009299.1(NM_017668.3):c.[41A>C;250del]",
        "normalized": "NG_009299.1(NM_017668.3):c.[41A>C;*189del]",
        "genomic": "NG_009299.1:g.[137591del;137800T>G]",
        "coding_protein_descriptions": {
            (
                "NG_009299.1(NM_001040113.2):c.[4138-31del;4316T>G]",
                "NG_009299.1(NP_001035202.1):p.(Leu1439Arg)",
                # TODO: Check if it should be '?'
            ),
            (
                "NG_009299.1(NM_002474.3):c.[4117-31del;4295T>G]",
                "NG_009299.1(NP_002465.1):p.(Leu1432Arg)",
                # TODO: Check if it should be '?'
            ),
            (
                "NG_009299.1(NM_001143979.2):c.[41A>C;*189del]",
                "NG_009299.1(NP_001137451.1):p.?",
                # TODO: Check if it should not be '?'
            ),
        },
        "protein_description": "NG_009299.1(NP_060138.1):p.?",
        # TODO: Check if it should not be '?'
        "to_test": True,
    },
    {
        "keywords": ["reverse strand"],
        "input": "NG_009299.1(NM_017668.3):c.[250del;41A>C]",
        "normalized": "NG_009299.1(NM_017668.3):c.[41A>C;*189del]",
        "infos": ["ICORRECTEDPOINT", "ISORTEDVARIANTS"],
        "to_test": True,
    },
    {
        "keywords": ["reverse strand"],
        "input": "NG_009299.1(NM_002474.3):c.[310del;295G>A]",
        "normalized": "NG_009299.1(NM_002474.3):c.[295G>A;311del]",
        "infos": ["ISORTEDVARIANTS"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.[274+20C>T;400_401insNM_003002.4:100_102]",
        "normalized": "NG_012337.1(NM_003002.2):c.[294C>T;399_401T[6]]",
        "infos": ["ICORRECTEDCOORDINATESYSTEM", "ICORRECTEDPOINT"],
        "to_test": True,
    },
    {
        "keywords": ["there should be no deleted sequence mismatch error"],
        "input": "NG_012337.1(NM_003002.2):c.2740000T>T",
        "infos": ["ICORRECTEDPOINT"],
        "errors": ["EOUTOFBOUNDARY"],
        "to_test": True,
    },
    {
        "keywords": ["no operation (equal)"],
        "input": "NM_002001.2:c.=",
        "normalized": "NM_002001.2:c.=",
        "infos": ["IWHOLETRANSCRIPTEXON"],
        "to_test": True,
    },
    {
        "keywords": ["no operation with selector (equal)"],
        "input": "NG_012337.1(NM_003002.2):c.=",
        "normalized": "NG_012337.1(NM_003002.2):c.=",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NM_003002.2:c.274+20C>T",
        "errors": ["EINTRONIC"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NM_003002.2:c.[100del;274+20C>T]",
        "errors": ["EINTRONIC"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NM_003002.2:c.[274+20C>T;300del]",
        "errors": ["EINTRONIC"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NM_003002.2:c.[100del;200_201insNM_003002.2:274+20]",
        "infos": [
            "ICORRECTEDCOORDINATESYSTEM",
            "IWHOLETRANSCRIPTEXON",
            "IWHOLETRANSCRIPTEXON",
        ],
        "errors": ["EINTRONIC"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NR_038420.1:n.100+10del",
        "errors": ["EINTRONIC"],
        "to_test": True,
    },
    {
        "keywords": ["ensembl", "plus strand"],
        "input": "ENST00000375549:c.100del",
        "normalized": "ENST00000375549.8:c.102del",
        "to_test": True,
    },
    {
        "keywords": ["ensembl", "plus strand"],
        "input": "ENST00000375549.8:c.100del",
        "normalized": "ENST00000375549.8:c.102del",
        "to_test": True,
    },
    {
        "keywords": ["ensembl", "plus strand"],
        "input": "ENSG00000204370.13(ENST00000375549.8):c.100del",
        "normalized": "ENSG00000204370.13(ENST00000375549.8):c.102del",
        "to_test": True,
    },
    {
        "keywords": ["ensembl"],
        "input": "ENST00000375549.7:c.100del",
        "errors": ["ERETR"],
        "to_test": True,
    },
    {
        "keywords": ["ensembl", "minus strand"],
        "input": "ENST00000452863.10:c.9del",
        "normalized": "ENST00000452863.10:c.10del",
        "to_test": True,
    },
    {
        "keywords": ["ensembl", "minus strand"],
        "input": "ENSG00000184937.16(ENST00000452863.10):c.9del",
        "normalized": "ENSG00000184937.16(ENST00000452863.10):c.10del",
        "to_test": True,
    },
    {
        "keywords": ["#46"],
        "input": "NG_012337.1:g.[104_105insA;105del;105_106insC]",
        "normalized": "NG_012337.1:g.105delinsAC",
        "to_test": True,
    },
    {
        "keywords": ["#46"],
        "input": "NG_012337.1:g.[104_105insA;105_106insC;105del]",
        "normalized": "NG_012337.1:g.105delinsAC",
        "infos": ["ISORTEDVARIANTS"],
        "to_test": True,
    },
    {
        "keywords": ["#46"],
        "input": "NG_012337.1:g.[105del;104_105insA;105_106insC]",
        "normalized": "NG_012337.1:g.105delinsAC",
        "infos": ["ISORTEDVARIANTS"],
        "to_test": True,
    },
    {
        "keywords": ["#46"],
        "input": "NG_012337.1:g.[105del;105_106insC;104_105insA]",
        "normalized": "NG_012337.1:g.105delinsAC",
        "infos": ["ISORTEDVARIANTS"],
        "to_test": True,
    },
    {
        "keywords": ["#46"],
        "input": "NG_012337.1:g.[105_106insC;104_105insA;105del]",
        "normalized": "NG_012337.1:g.105delinsAC",
        "infos": ["ISORTEDVARIANTS"],
        "to_test": True,
    },
    {
        "keywords": ["#46"],
        "input": "NG_012337.1:g.[105_106insC;105del;104_105insA]",
        "normalized": "NG_012337.1:g.105delinsAC",
        "infos": ["ISORTEDVARIANTS"],
        "to_test": True,
    },
    {
        "keywords": ["#46", "reverse strand"],
        "input": "NG_012337.1(NM_012459.2):c.[7_8insC;8del]",
        "normalized": "NG_012337.1(NM_012459.2):c.8A>C",
        "to_test": True,
    },
    {
        "keywords": ["#46", "reverse strand"],
        "input": "NG_012337.1(NM_012459.2):c.[8del;7_8insC]",
        "normalized": "NG_012337.1(NM_012459.2):c.8A>C",
        "infos": ["ISORTEDVARIANTS"],
        "to_test": True,
    },
    {
        "keywords": ["#49", "mitochondrion"],
        "input": "NC_012920.1:g.3243A>G",
        "normalized": "NC_012920.1:m.3243A>G",
        "infos": ["ICORRECTEDCOORDINATESYSTEM"],
        "to_test": True,
    },
    {
        "keywords": ["#49", "mitochondrion"],
        "input": "NC_012920.1:m.3243A>G",
        "normalized": "NC_012920.1:m.3243A>G",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.[4_5insT;4insA;4_5insA]",
        "errors": ["EINSERTIONRANGE", "EOVERLAP"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.[4_5insT;4_5insA]",
        "errors": ["EOVERLAP"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.[4del;4_5insT;4insA;4_5insA]",
        "errors": ["EINSERTIONRANGE", "EOVERLAP"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.[4_5insT;4insA;4_5insA;5del]",
        "errors": ["EINSERTIONRANGE", "EOVERLAP"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.[4del;4_5insT;4insA;4_5insA;5del]",
        "errors": ["EINSERTIONRANGE", "EOVERLAP"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.[273_274insT;274G>T;274_275insA]",
        "normalized": "NG_012337.1(NM_003002.2):c.274delinsTTA",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.[274_275delinsT;274_275insA]",
        "errors": ["EOVERLAP"],
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.[273_274insT;274G>T;274_275insA]",
        "normalized": "NG_012337.1(NM_003002.2):c.274delinsTTA",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.pterdel",
        "normalized": "NG_012337.1(NM_003002.2):c.-5059del",
        "genomic": "NG_012337.1:g.3del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.pterdel",
        "normalized": "NG_012337.1:g.3del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.qterdel",
        "normalized": "NG_012337.1(NM_003002.2):c.*2824del",
        "genomic": "NG_012337.1:g.15948del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NM_003002.2:c.pterdel",
        "normalized": "NM_003002.2:c.-61del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NM_003002.2:c.qterdel",
        "normalized": "NM_003002.2:c.*841del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "LRG_24:g.5525_5532delinsNM_003002.2:c.pter_-51",
        "normalized": "LRG_24:g.5525_5532delinsGTGGGAATTGT",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "LRG_24:g.5525_5532delinsNM_003002.2:c.*835_qter",
        "normalized": "LRG_24:g.5525_5532delinsAAAAAAA",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "LRG_24(t1):c.126_133delins[NM_003002.2:c.pter_-51;NM_003002.2:c.*835_qter]",
        "normalized": "LRG_24(t1):c.126_133delinsGTGGGAATTGTAAAAAAA",
        "genomic": "LRG_24:g.5525_5532delinsGTGGGAATTGTAAAAAAA",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "LRG_24(t1):c.pter_qterdelinspter_qter",
        "normalized": "LRG_24(t1):c.=",
        "genomic": "LRG_24:g.=",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "LRG_24(t1):c.pter_qterdelins[pter_qter;NM_003002.2:c.*835_qter]",
        "normalized": "LRG_24(t1):c.*2327_*2328insAAAAAAA",
        "genomic": "LRG_24:g.11486_11487insAAAAAAA",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.[274=;275del]",
        "normalized": "NG_012337.1(NM_003002.2):c.275del",
        "genomic": "NG_012337.1:g.7126del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_003002.2):c.[274del;275=]",
        "normalized": "NG_012337.1(NM_003002.2):c.274del",
        "genomic": "NG_012337.1:g.7125del",
        "to_test": True,
    },
    {
        "keywords": ["repeat", "dbsnp"],
        "input": "NG_012337.1:g.4917GC[5]",
        "genomic": "NG_012337.1:g.4917_4920GC[5]",
        "to_test": True,
    },
    {
        "keywords": ["repeat", "dbsnp", "reverse"],
        "input": "NG_012337.1(NM_012459.2):c.3GC[5]",
        "normalized": "NG_012337.1(NM_012459.2):c.3_6GC[5]",
        "protein_description": "NG_012337.1(NP_036591.2):p.(Arg2_Lys3insAlaArg)",
        "genomic": "NG_012337.1:g.4917_4920GC[5]",
        "to_test": True,
    },
    {
        "keywords": ["repeat", "dbsnp"],
        "input": "NG_012337.1:g.4911GT[5]",
        "genomic": "NG_012337.1:g.4911_4914GT[5]",
        "to_test": True,
    },
    {
        "keywords": ["repeat", "dbsnp", "reverse"],
        "input": "NG_012337.1(NM_012459.2):c.10CA[5]",
        "normalized": "NG_012337.1(NM_012459.2):c.10_13CA[5]",
        "genomic": "NG_012337.1:g.4911_4914GT[5]",
        "to_test": True,
    },
    {
        "keywords": ["duplication", "reverse"],
        "input": "NG_012337.1(NM_012459.2):c.5_6dup",
        "normalized": "NG_012337.1(NM_012459.2):c.5_6dup",
        "genomic": "NG_012337.1:g.4919_4920dup",
        "to_test": True,
    },
    {
        "keywords": ["duplication", "reverse"],
        "input": "NG_012337.1(NM_012459.2):c.5dup",
        "normalized": "NG_012337.1(NM_012459.2):c.5dup",
        "genomic": "NG_012337.1:g.4918dup",
        "to_test": True,
    },
    {
        "keywords": ["superfluous location should be within sequence (#69)"],
        "input": "NM_003002.4:c.[-100A>C;20000A>C]",
        "errors": [
            "EOUTOFBOUNDARY",
            "EOUTOFBOUNDARY",
        ],
        "to_test": True,
    },
    {
        "keywords": ["protein reverse strand exon end"],
        "input": "NG_008835.1(NM_022153.2):c.568del",
        "normalized": "NG_008835.1(NM_022153.2):c.568del",
        "genomic": "NG_008835.1:g.368924del",
        "coding_protein_descriptions": {
            (
                "NG_008835.1(NM_022124.6):c.4846-16810del",
                "NG_008835.1(NP_071407.4):p.(=)",
            ),
        },
        "protein_description": "NG_008835.1(NP_071436.1):p.(Asn190Thrfs*14)",
        "rna_description": "NG_008835.1(NM_022153.2):r.(569del)",
        "to_test": True,
    },
    {
        "keywords": ["protein reverse strand exon end splice site error"],
        "input": "NG_008835.1(NM_022153.2):c.82del",
        "normalized": "NG_008835.1(NM_022153.2):c.82+1del",
        "genomic": "NG_008835.1:g.381412del",
        "coding_protein_descriptions": {
            (
                "NG_008835.1(NM_022124.6):c.4846-4322del",
                "NG_008835.1(NP_071407.4):p.(=)",
            ),
        },
        "to_test": True,
    },
    # {
    #     "keywords": [
    #
    #     ],
    #     "input": "",
    #     "normalized": "",
    #     "genomic": "",
    #     "to_test": False,
    # },
]

TESTS_ALL = M2_TESTS + TESTS


OTHER = [
    ("NG_012337.1(SDHD_v001):c.274G>T", "NG_012337.1:g.7125G>T"),
    # ('NG_012337.1(NM_003002.2):c.274G>T',
    #  'NG_012337.1:g.7125G>T'),
    # Roll
    # ----
    # We should roll.
    ("NM_003002.4:c.273del", "NM_003002.4:n.309del"),
    # We cannot roll, we are at the end.
    ("NM_003002.4:c.274del", "NM_003002.4:n.309del"),
    # We can roll but should not, because it is over a splice site.
    ("NM_000088.3:c.333del", "NM_000088.3:n.459del"),
    # We can roll two positions, but should roll only one because
    # otherwise it is over a splice site.
    ("NM_000088.3:c.368del", "NM_000088.3:n.495del"),
    # We can roll and should, we stay in the same exon.
    ("NM_000088.3:c.334del", "NM_000088.3:n.461del"),
    # ('NG_012337.1:c.274G>T',
    #  'NG_012337.1:g.7125G>T'),
    # ('NC_000024.10(PRY_v001):c.948delC',
    #  'NC_000024.10:g.22514574del'),
    # ('NC_000024.10(XR_001756027.1):c.948del',
    #  'NC_000024.10:g.1352647del'),
    # ('NC_000001.11(ZRANB2-AS1_v001):c.948delC',
    #  'NC_000001.11:g.71047454del'),
    # ('NC_000001.11(NR_038420.1):c.948delC',
    #  'NC_000001.11:g.71047454del'),
    # ('NG_009497.1(LOC102723833_v001):c.8682-19dup',
    #  'NG_009497.1:g.561208dup'),
    # ('NM_003002.4:c.5762_5763insNG_009113.2:g.91138_91274',
    #  'NG_009497.1:g.561208dup'),
    # ('NG_009113.2:g.575_576insNG_009113.2:g.9118_9127',
    #  'NG_009113.2:g.575_576insGGAGGCAGAG'),
    # ('NG_009113.2:g.575_576insNG_009497.1:g.9118_9127',
    #  'NG_009113.2:g.575_576insAGTTGGACTG'),
    # ('NG_017013.2(NM_001126118.1):c.100-10del',
    #  'Offset may not be from position 100 because this is not an exon boundary.'),
    # ('NG_017013.2(NM_001126118.1):c.259del',
    #  ''),
]
