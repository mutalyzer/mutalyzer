

TESTS_ALL = [
    {
        'keywords': ['correct input', 'no errors', 'no warnings',
                     'genomic reference', 'NCBI', 'g.',
                     'single variant', 'substitution'],
        'input': 'NG_012337.1:g.4C>T',
        'normalized': 'NG_012337.1:g.4C>T',
        'genomic': 'NG_012337.1:g.4C>T',
        'to_test': True
    },
    {
        'keywords': ['correct input', 'no errors', 'no warnings',
                     'genomic reference', 'NCBI', 'g.',
                     'single variant', 'deletion range'],
        'input': 'NG_017013.2:g.17013_17014del',
        'normalized': 'NG_017013.2:g.17013_17014del',
        'to_test': True
    },
    {
        'keywords': ['correct input', 'no errors', 'no warnings',
                     'genomic reference', 'NCBI', 'g.',
                     'single variant', 'deletion insertion single range'],
        'input': 'NG_012337.1:g.4delins7_31',
        'normalized': 'NG_012337.1:g.4delins7_31',
        'to_test': True
    },
    {
        'keywords': ['correct input', 'no errors', 'no warnings',
                     'genomic reference', 'NCBI', 'g.',
                     'single variant', 'duplication single'],
        'input': 'NG_017013.2:g.19258dup',
        'normalized': 'NG_017013.2:g.19258dup',
        'to_test': True
    },
    {
        'keywords': ['correct input', 'no errors', 'no warnings',
                     'genomic reference', 'NCBI', 'g.',
                     'single variant', 'duplication range'],
        'input': 'NG_017013.2:g.16508_16509dup',
        'normalized': 'NG_017013.2:g.16508_16509dup',
        'to_test': True
    },
    {
        'keywords': ['correct input', 'no errors', 'no warnings',
                     'genomic reference', 'NCBI', 'g.',
                     'single variant', 'insertion single sequence'],
        'input': 'NG_012337.1:g.90_91insC',
        'normalized': 'NG_012337.1:g.90_91insC',
        'to_test': True
    },
    # ---------
    {
        'keywords': ['correct input', 'no errors', 'no warnings',
                     'NCBI', 'c.', 'single variant', 'substitution',
                     'roll forward'],
        'input': 'NG_007485.1(NM_000077.4):c.274G>T',
        'normalized': 'NG_007485.1(NM_000077.4):c.274G>T',
        'genomic': 'NG_007485.1:g.28407G>T',
        'to_test': True
    },
    # ---------
    {
        'keywords': ['corrected', 'no errors', 'warnings',
                     'NCBI', 'g.', 'single variant', 'deletion range',
                     'roll forward'],
        'input': 'NG_012337.1:g.26_31del',
        'normalized': 'NG_012337.1:g.29_34del',
        'warnings': 'Deletion at position 26_31 was given, however, the HGVS '
                    'notation prescribes that on the forward strand it should '
                    'be at position 29_34.',
        'to_test': True
    },
    {
        'keywords': ['corrected', 'no errors', 'warnings',
                     'NCBI', 'g.', 'single variant', 'deletion range',
                     'roll forward'],
        'input': 'NG_017013.2:g.17011_17012del',
        'normalized': 'NG_017013.2:g.17013_17014del',
        'to_test': True
    },
    # ---------

    {
        'keywords': [],
        'input': 'NG_017013.2:g.17415_17417delinsGCG',
        'normalized': 'NG_017013.2:g.[17415C>G;17417A>G]',
        'genomic': 'NG_017013.2:g.[17415C>G;17417A>G]',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_017013.2:g.17496_17497insAGCTGCTCAGATAGCGA',
        'normalized': 'NG_017013.2:g.17496_17497ins[A;17481_17496]',
        'genomic': '',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_017013.2:g.1_32772del',
        'normalized': 'NG_017013.2:g.1_32772del',
        'genomic': '',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_017013.2:g.[16985A>T;17013_17014del]',
        'normalized': 'NG_017013.2:g.[16985A>T;17013_17014del]',
        'genomic': '',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_017013.2:g.[16985A>T;17011_17012del]',
        'normalized': 'NG_017013.2:g.[16985A>T;17013_17014del]',
        'genomic': '',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_017013.2:g.17011_17012del',
        'normalized': 'NG_017013.2:g.17013_17014del',
        'genomic': '',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_012337.1:g.4delins7_50',
        'normalized': 'NG_012337.1:g.[3_4insGGTT;5_6ins12_50]',
        'genomic': 'NG_012337.1:g.[3_4insGGTT;5_6ins12_50]',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_012337.1:g.26_31del',
        'normalized': 'NG_012337.1:g.29_34del',
        'genomic': '',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_012337.1:g.100_200delins100_101',
        'normalized': 'NG_012337.1:g.102_200del',
        'genomic': 'NG_012337.1:g.102_200del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_017013.2:g.17471_17471del',
        'normalized': 'NG_017013.2:g.17471del',
        'genomic': '',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_017013.2:g.18748_18750delinsCAT',
        'normalized': 'NG_017013.2:g.18749G>A',
        'genomic': '',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_017013.2:g.17394_17395insC',
        'normalized': 'NG_017013.2:g.17394dup',
        'genomic': 'NG_017013.2:g.17394dup',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_012337.1(NM_003002.2):c.274G>T',
        'normalized': 'NG_012337.1(NM_003002.2):c.274G>T',
        'genomic': 'NG_012337.1:g.7125G>T',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_012337.1(NM_003002.2):c.1del',
        'normalized': 'NG_012337.1(NM_003002.2):c.1del',
        'genomic': 'NG_012337.1:g.5062del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_012337.1(NM_003002.2):c.-1del',
        'normalized': 'NG_012337.1(NM_003002.2):c.-1del',
        'genomic': 'NG_012337.1:g.5061del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_012337.1(NM_003002.2):c.52+1del',
        'normalized': 'NG_012337.1(NM_003002.2):c.52+1del',
        'genomic': 'NG_012337.1:g.5114del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_012337.1(NM_003002.2):c.*824del',
        'normalized': 'NG_012337.1(NM_003002.2):c.*824del',
        'genomic': 'NG_012337.1:g.13948del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_012337.1(NM_003002.2):c.*824+10del',
        'normalized': 'NG_012337.1(NM_003002.2):c.*834del',
        'genomic': 'NG_012337.1:g.13958del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NM_003002.4:c.1del',
        'normalized': 'NM_003002.4:c.1del',
        'coordinate': 'NM_003002.4:x.36del',
        'to_test': False
    },
    {
        'keywords': [],
        'input': 'NG_029724.1(NM_004321.7):c.101del',
        'normalized': 'NG_029724.1(NM_004321.7):c.102del',
        'genomic': 'NG_029724.1:g.27557del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_029724.1(NM_004321.7):c.101del',
        'normalized': '',
        'genomic': 'NG_029724.1:g.27557del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_029724.1:g.10_20del',
        'normalized': 'NG_029724.1:g.10_20del',
        'genomic': 'NG_029724.1:g.10_20del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_029724.1:g.10_20del11',
        'normalized': 'NG_029724.1:g.10_20del',
        'genomic': 'NG_029724.1:g.10_20del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_029724.1:g.10delG',
        'normalized': 'NG_029724.1:g.10del',
        'genomic': 'NG_029724.1:g.10del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_029724.1:g.10_11delGT',
        'normalized': 'NG_029724.1:g.10_11del',
        'genomic': 'NG_029724.1:g.10_11del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_008835.1(CDH23_v001):c.1449+846delA',
        'normalized': '',
        'genomic': 'NG_008835.1:g.255529del',
        'to_test': False
    },
    {
        'keywords': [],
        'input': 'NG_009113.2(NR2E3_v001):c.948delC',
        'normalized': '',
        'genomic': 'NG_009113.2:g.8038del',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_007107.2(MECP2_v001):c.378-17delT',
        'normalized': '',
        'genomic': 'NG_007107.2:g.110661del',
        'to_test': False
    },
    {
        'keywords': [],
        'input': 'NG_009113.2(NR2E3_v002):c.948delC',
        'normalized': '',
        'genomic': 'NG_009113.2:g.8038del',
        'to_test': False
    },
    {
        'keywords': [],
        'input': 'NG_009113.2(NM_016346.4):c.948delC',
        'normalized': '',
        'genomic': 'NG_009113.2:g.8038del',
        'to_test': False
    },
    {
        'keywords': [],
        'input': 'NG_009497.1(NM_206933.2):c.8682-19dupT',
        'normalized': 'NG_009497.1(NM_206933.2):c.8682-19dup',
        'genomic': 'NG_009497.1:g.561208dup',
        'to_test': True
    },
    {
        'keywords': [],
        'input': 'NG_009497.1(USH2A_v001):c.8682-19dup',
        'normalized': '',
        'genomic': 'NG_009497.1:g.561208dup',
        'to_test': False
    },
    # {
    #     'keywords': [],
    #     'input': '',
    #     'normalized': '',
    #     'genomic': '',
    #     'to_test': False
    # },
]