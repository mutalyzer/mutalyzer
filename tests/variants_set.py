from .m2 import TESTS

TESTS_ALL = [
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
        "genomic": "NG_012337.1:g.4C>T",
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
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.1_32772del",
        "normalized": "NG_017013.2:g.1_32772del",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.[16985A>T;17013_17014del]",
        "normalized": "NG_017013.2:g.[16985A>T;17013_17014del]",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.[16985A>T;17011_17012del]",
        "normalized": "NG_017013.2:g.[16985A>T;17013_17014del]",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.17011_17012del",
        "normalized": "NG_017013.2:g.17013_17014del",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.4delins7_50",
        "normalized": "NG_012337.1:g.[3_4insGGTT;5_6insACCATATCTCTACTTTGTGTTTATGTTTGTGTATGCATT]",
        "genomic": "NG_012337.1:g.[3_4insGGTT;5_6insACCATATCTCTACTTTGTGTTTATGTTTGTGTATGCATT]",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1:g.26_31del",
        "normalized": "NG_012337.1:g.29_34del",
        "genomic": "",
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
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_017013.2:g.18748_18750delinsCAT",
        "normalized": "NG_017013.2:g.18749G>A",
        "genomic": "",
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
        "coordinate": "NM_003002.4:x.36del",
        "to_test": False,
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
        "genomic": "NG_029724.1:g.10_20del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_029724.1:g.10_20del11",
        "normalized": "NG_029724.1:g.10_20del",
        "genomic": "NG_029724.1:g.10_20del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_029724.1:g.10delG",
        "normalized": "NG_029724.1:g.10del",
        "genomic": "NG_029724.1:g.10del",
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_029724.1:g.10_11delGT",
        "normalized": "NG_029724.1:g.10_11del",
        "genomic": "NG_029724.1:g.10_11del",
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
        "normalized": "",
        "genomic": "NG_008835.1:g.255529del",
        "to_test": False,
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
        "input": "NG_007107.2(MECP2_v001):c.378-17delT",
        "normalized": "",
        "genomic": "NG_007107.2:g.110661del",
        "to_test": False,
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
        "normalized": "",
        "genomic": "NG_009113.2:g.8038del",
        "to_test": False,
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
        "normalized": "",
        "genomic": "NG_009497.1:g.561208dup",
        "to_test": False,
    },
    {
        "keywords": [],
        "input": "NG_007485.1(NM_058195.3):c.141_142del",
        "normalized": "NG_007485.1(NM_058195.3):c.141_142del",
        "genomic": "NG_007485.1:g.5301_5302del",
        "coding_protein_descriptions": {
            (
                "NG_007485.1(NM_058195.3):c.141_142del",
                "NG_007485.1(NP_478102.2):p.(Met48Alafs*14)",
            ),
            (
                "NG_007485.1(NM_000077.4):c.-19364_-19363del",
                "NG_007485.1(NP_000068.1):p.(=)",
            ),
        },
        "noncoding": ["NG_007485.1(NR_024274.1):n.616+26436_616+26437del"],
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6932_6933insAGCAACGTGATCGCCTCCCTCACCT",
        "normalized": "LRG_303:g.6908_6932dup",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6932_6933insAGCAACGTGATCGCCTCCCTCACCT[0]",
        "normalized": "LRG_303:g.1_11312=",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6932_6933insAGCAACGTGATCGCCTCCCTCACCT[1]",
        "normalized": "LRG_303:g.6908_6932dup",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6932_6933insAGCAACGTGATCGCCTCCCTCACCT[2]",
        "normalized": "LRG_303:g.6908_6932AGCAACGTGATCGCCTCCCTCACCT[3]",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.6908_6932AGCAACGTGATCGCCTCCCTCACCT[0]",
        "normalized": "LRG_303:g.6908_6932del",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.4_5insGT[5]",
        "normalized": "LRG_303:g.1_4GT[7]",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.3_4GT[5]",
        "normalized": "LRG_303:g.1_4GT[6]",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.1_2GT[0]",
        "normalized": "LRG_303:g.3_4del",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.3_4GT[0]",
        "normalized": "LRG_303:g.3_4del",
        "genomic": "",
        "to_test": True,
    },
    {
        "keywords": ["repeat"],
        "input": "LRG_303:g.3_4GT[2]",
        "normalized": "LRG_303:g.3_4dup",
        "normalized-alt": "LRG_303:g.1_4[3]",
        "genomic": "",
        "to_test": False,
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
                "NG_009299.1(NM_017668.3):c.41A>C",
                "NG_009299.1(NP_060138.1):p.(Gln14Profs*68)",
            ),
            (
                "NG_009299.1(NM_001143979.2):c.41A>C",
                "NG_009299.1(NP_001137451.1):p.(Gln14Profs*68)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "",
        "normalized": "",
        "genomic": "",
        "to_test": False,
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
        "keywords": ["protein", "M2: similar to `accno_as_transcript_variant`"],
        "input": "NG_012337.1(NM_012459.2):c.4_5insGTA",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_012459.2):c.5_6insTAG",
                "NG_012337.1(NP_036591.2):p.(Arg2_Lys3insSer)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_012459.2):c.5_6delinsTAG",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_012459.2):c.5_6delinsTAG",
                "NG_012337.1(NP_036591.2):p.(Arg2Leufs*23)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": [],
        "input": "NG_012337.1(NM_012459.2):c.4_6delinsGTA",
        "coding_protein_descriptions": {
            (
                "NG_012337.1(NM_012459.2):c.4_6delinsGTA",
                "NG_012337.1(NP_036591.2):p.(Arg2Val)",
            ),
        },
        "to_test": True,
    },
    {
        "keywords": ["M2: `accno_as_transcript_variant` the error part"],
        "input": "NG_012337.1(DUMMYACCNO_9999.9):c.12_13insGATC",
        "errors": ["ENOSELECTORFOUND"],
        "to_test": True,
    },
    {
        "keywords": [
            "M2: dup_range_arg_seq_reverse(",
            "Range duplication with sequence argument (reverse strand).",
            "Switched from AL449423.14 to NG_007485.1 and from CDKN2A_v001 to NM_000077.4",
        ],
        # "input": "AL449423.14(CDKN2A_v001):c.41_44dupACTG",
        "input": "NG_007485.1(NM_000077.4):c.41_44dupACTG",
        "normalized": "NG_007485.1(NM_000077.4):c.41_44dup",
        "genomic": "NG_007485.1:g.24705_24708dup",
        # TODO: Check the protein descriptions.
        # TODO: Add test for dup sequence.
        "to_test": True,
    },
    {
        "keywords": ["reference", "LRG", "replace"],
        "input": "LRG_303(t1):c.10_11insLRG_1t1:c.100_101",
        "normalized": "LRG_303(t1):c.10_11insGA",
        "infos": ["ICORRECTEDLRGREFERENCE"],
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
] + TESTS


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
    #  'NG_009497.1:g.561208dup'),
    # ('NG_009113.2:g.575_576insNG_009497.1:g.9118_9127',
    #  'NG_009497.1:g.561208dup'),
    # ('NG_017013.2(NM_001126118.1):c.100-10del',
    #  'Offset may not be from position 100 because this is not an exon boundary.'),
    # ('NG_017013.2(NM_001126118.1):c.259del',
    #  'Offset may not be from position 100 because this is not an exon boundary.'),
]
