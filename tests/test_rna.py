import pytest
from mutalyzer_hgvs_parser import to_model

from mutalyzer.checker import splice_sites
from mutalyzer.converter.to_rna import (
    _trim_to_exons,
    get_location_type,
    get_position_type,
    to_rna_reference_model,
    to_rna_variants,
)
from mutalyzer.normalizer import normalize
from mutalyzer.reference import get_internal_selector_model, retrieve_reference

from .commons import code_in, patch_retriever

ESPLICESITE = {
    "code": "ESPLICESITE",
    "details": "Splice site(s) affected.",
    "paths": [["variants", 0]],
}

IVARIANTDISCARDED = {
    "code": "IVARIANTDISCARDED",
    "details": "Variant discarded.",
    "paths": [["variants", 0]],
}

CORRECTEDSEQUENCE = {
    "code": "CORRECTEDSEQUENCE",
    "details": "CORRECTEDSEQUENCE",
    "details": 'Sequence "AAA" corrected to "aaa".',
}

ISPLICESITEREMOVED = {
    "code": "ISPLICESITEREMOVED",
    "details": "Splice(s) sites removed.",
    "paths": [["variants", 0]],
}

IWHOLETRANSCRIPTEXON = {
    "code": "IWHOLETRANSCRIPTEXON",
    "details": "No exon features found in the 'NM_003002.2' reference sequence "
    "for 'NM_003002.2'. The entire transcript was assumed as one "
    "exon.",
    "paths": [("reference", "selector", "id")],
}

TESTS = [
    {
        "keywords": ["rna", "equal", "genomic", "mRNA"],
        "input": "NG_012337.1(NM_003002.2):r.275=",
        "normalized": "NG_012337.1(NM_003002.2):r.275=",
        "to_test": True,
    },
    {
        "keywords": ["rna", "uncertain", "genomic", "mRNA"],
        "input": "NG_012337.1(NM_003002.2):r.?",
        "normalized": "NG_012337.1(NM_003002.2):r.?",
        "to_test": False,
    },
    {
        "keywords": ["rna", "substitution", "genomic", "mRNA"],
        "input": "NG_012337.1(NM_003002.2):r.275a>c",
        "normalized": "NG_012337.1(NM_003002.2):r.275a>c",
        "protein_description": "NG_012337.1(NP_002993.1):p.(Asp92Ala)",
        "to_test": True,
    },
    {
        "keywords": ["rna", "substitution", "mRNA"],
        "input": "NM_003002.2:r.275a>c",
        "normalized": "NM_003002.2:r.275a>c",
        "protein_description": "NM_003002.2(NP_002993.1):p.(Asp92Ala)",
        "infos": [IWHOLETRANSCRIPTEXON],
        "to_test": True,
    },
    {
        "keywords": ["rna", "substitution", "ncRNA"],
        "input": "NR_038420.1:r.75g>a",
        "normalized": "NR_038420.1:r.75g>a",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "genomic", "mRNA"],
        "input": "NG_012337.1(NM_003002.2):r.273del",
        "normalized": "NG_012337.1(NM_003002.2):r.274del",
        "protein_description": "NG_012337.1(NP_002993.1):p.(Asp92Thrfs*43)",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "genomic", "mRNA", "same exon", "plus strand"],
        "input": "NG_012337.3(NM_003002.4):r.-30_40del",
        "normalized": "NG_012337.3(NM_003002.4):r.-30_40del",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "mRNA", "same exon"],
        "input": "NM_003002.4:r.-30_40del",
        "normalized": "NM_003002.4:r.-30_40del",
        "protein_description": "NM_003002.4(NP_002993.1):p.?",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "genomic", "mRNA", "same exon", "minus strand"],
        "input": "NG_012337.3(NM_012459.4):r.-30_40del",
        "normalized": "NG_012337.3(NM_012459.4):r.-30_40del",
        "protein_description": "NG_012337.3(NP_036591.3):p.?",
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "mRNA", "same exon"],
        "input": "NM_012459.4:r.-30_40del",
        "normalized": "NM_012459.4:r.-30_40del",
        "protein_description": "NM_012459.4(NP_036591.3):p.?",
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "same exon",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_40delinsAAA",
        "normalized": "NG_012337.3(NM_003002.4):r.-30_40delinsaaa",
        "protein_description": "NG_012337.3(NP_002993.1):p.?",
        "infos": [CORRECTEDSEQUENCE],
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion insertion", "mRNA", "same exon"],
        "input": "NM_003002.4:r.-30_40delinsAAA",
        "normalized": "NM_003002.4:r.-30_40delinsaaa",
        "protein_description": "NM_003002.4(NP_002993.1):p.?",
        "infos": [CORRECTEDSEQUENCE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "same exon",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_40delinsAAA",
        "normalized": "NG_012337.3(NM_012459.4):r.-30_40delinsaaa",
        "protein_description": "NG_012337.3(NP_036591.3):p.?",
        "infos": [CORRECTEDSEQUENCE],
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion insertion", "mRNA", "same exon"],
        "input": "NM_012459.4:r.-30_40delinsAAA",
        "normalized": "NM_012459.4:r.-30_40delinsaaa",
        "protein_description": "NM_012459.4(NP_036591.3):p.?",
        "infos": [CORRECTEDSEQUENCE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - exon",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_200del",
        "normalized": "NG_012337.3(NM_003002.4):r.-30_200del",
        "protein_description": "NG_012337.3(NP_002993.1):p.?",
        "infos": [ISPLICESITEREMOVED],
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "mRNA", "exon - exon"],
        "input": "NM_003002.4:r.-30_200del",
        "normalized": "NM_003002.4:r.-30_200del",
        "protein_description": "NM_003002.4(NP_002993.1):p.?",
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - exon",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_200del",
        "normalized": "NG_012337.3(NM_012459.4):r.-29_201del",
        "protein_description": "NG_012337.3(NP_036591.3):p.?",
        "infos": [ISPLICESITEREMOVED],
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion", "mRNA", "exon - exon"],
        "input": "NM_012459.4:r.-30_200del",
        "normalized": "NM_012459.4:r.-29_201del",
        "protein_description": "NM_012459.4(NP_036591.3):p.?",
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - exon",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_200delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion insertion", "mRNA", "exon - exon"],
        "input": "NM_003002.4:r.-30_200delinsAAA",
        "normalized": "NM_003002.4:r.-30_200delinsaaa",
        "protein_description": "NM_003002.4(NP_002993.1):p.?",
        "infos": [CORRECTEDSEQUENCE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - exon",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_200delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": ["rna", "deletion insertion", "mRNA", "exon - exon"],
        "input": "NM_012459.4:r.-30_200delinsAAA",
        "normalized": "NM_012459.4:r.-30_200delinsaaa",
        "protein_description": "NM_012459.4(NP_036591.3):p.?",
        "infos": [CORRECTEDSEQUENCE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_52+1del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_53-1del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_169+1del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_84+1del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_85-1del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_52+2delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_53-2delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_169+2delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_84+2delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_85-2delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - around splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_52+3del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - around splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_53-3del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - around splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_169+3del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - around splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_84+3del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - around splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_85-3del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - around splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_52+7delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - around splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_53-7delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - around splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_84+7delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "exon - around splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_85-7delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - intron",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_52+8del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - intron",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_53-8del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - intron",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.-30_169+8del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "exon - intron",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.-30_84+8del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    # ---- reduced the number of checks
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+1_52+2del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+1_84+2del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+1_53-2del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+1_85-2del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - around splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+1_52+3del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - around splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+1_84+3del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - around splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+1_53-3del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - around splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+1_85-3del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - intron",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+1_52+8del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "splice site - intron",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+1_84+8del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "around splice site - around splice site",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+3_52+7del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "around splice site - around splice site",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+3_84+7del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "around splice site - intron",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+3_52+8del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "around splice site - intron",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+3_84+8del",
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "same intron",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+8_53-8del",
        "normalized": "NG_012337.3(NM_003002.4):r.=",
        "infos": [IVARIANTDISCARDED],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "intron - intron",
            "same intron",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+8_85-8del",
        "normalized": "NG_012337.3(NM_012459.4):r.=",
        "infos": [IVARIANTDISCARDED],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "intron - intron",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+8_169+8del",
        "normalized": "NG_012337.3(NM_003002.4):r.55_171del",
        "protein_description": "NG_012337.3(NP_002993.1):p.(Leu19_Ser57del)",
        "infos": [ISPLICESITEREMOVED],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion",
            "genomic",
            "mRNA",
            "intron - intron",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+8_*503del",
        "normalized": "NG_012337.3(NM_012459.4):r.85_*495del",
        "infos": [ISPLICESITEREMOVED],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "intron - intron",
            "plus strand",
        ],
        "input": "NG_012337.3(NM_003002.4):r.52+8_169+8delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "deletion insertion",
            "genomic",
            "mRNA",
            "intron - intron",
            "minus strand",
        ],
        "input": "NG_012337.3(NM_012459.4):r.84+8_*503delinsAAA",
        "infos": [CORRECTEDSEQUENCE],
        "errors": [ESPLICESITE],
        "to_test": True,
    },
    # ---- other
    {
        "keywords": [
            "rna",
            "substitution",
            "mRNA",
        ],
        "input": "NM_003002.2:r.277g>u",
        "infos": [IWHOLETRANSCRIPTEXON],
        "errors": [
            {
                "code": "ESEQUENCEMISMATCH",
                "details": "g not found in the reference sequence, found u instead.",
                "paths": [["variants", 0, "deleted"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "substitution",
            "mRNA",
        ],
        "input": "NM_003002.2:r.277t>u",
        "infos": [IWHOLETRANSCRIPTEXON],
        "errors": [
            {
                "code": "ENORNA",
                "details": "Sequence t is not an RNA sequence.",
                "paths": [["variants", 0, "deleted", 0, "sequence"]],
            }
        ],
        "to_test": True,
    },
    {
        "keywords": [
            "rna",
            "allele",
        ],
        "input": "NG_012337.1(NM_003002.2):r.([274g>u;278u>g])",
        "errors": [
            {
                "code": "ESEQUENCEMISMATCH",
                "details": "u not found in the reference sequence, found a instead.",
                "paths": [["variants", 1, "deleted"]],
            }
        ],
        "to_test": True,
    },
]


def get_tests(tests, t_type):
    output = []
    for test in tests:
        if test.get("to_test"):
            if test.get(t_type) is not None:
                output.append((test["input"], test[t_type]))
            else:
                output.append((test["input"], None))
    return output


@pytest.mark.parametrize(
    "input_description, normalized",
    get_tests(TESTS, "normalized"),
)
def test_rna(input_description, normalized):
    d = normalize(input_description)
    assert d.get("normalized_description") == normalized


@pytest.mark.parametrize(
    "input_description, errors",
    get_tests(TESTS, "errors"),
)
def test_rna_errors(input_description, errors):
    d = normalize(input_description)
    assert d.get("errors") == errors


@pytest.mark.parametrize(
    "input_description, infos",
    get_tests(TESTS, "infos"),
)
def test_rna_infos(input_description, infos):
    d = normalize(input_description)
    assert d.get("infos") == infos


@pytest.mark.parametrize(
    "input_description, protein_description",
    get_tests(TESTS, "protein_description"),
)
def test_rna_protein(input_description, protein_description):

    normalized_output = normalize(input_description)
    if protein_description is not None:
        normalizer_protein = normalized_output["protein"]["description"]
        assert normalizer_protein == protein_description


@pytest.mark.parametrize(
    "r_id, s_id, expected",
    [
        (
            "TEST_REF",
            "NM_PLUS",
            {
                "annotations": {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 0},
                        "end": {"type": "point", "position": 636},
                    },
                    "id": "TEST_REF",
                    "type": "record",
                    "qualifiers": {
                        "name": "TEST_REF",
                        "chromosome": "TEST_CHR",
                        "mol_type": "genomic DNA",
                    },
                    "features": [
                        {
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": 0},
                                "end": {"type": "point", "position": 636},
                                "strand": 1,
                            },
                            "id": "PLUS",
                            "type": "gene",
                            "qualifiers": {"name": "PLUS"},
                            "features": [
                                {
                                    "id": "NM_PLUS",
                                    "type": "mRNA",
                                    "location": {
                                        "type": "range",
                                        "start": {"type": "point", "position": 0},
                                        "end": {"type": "point", "position": 636},
                                        "strand": 1,
                                    },
                                    "features": [
                                        {
                                            "id": "exon-NM_PLUS-1",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 54,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                        {
                                            "id": "exon-NM_PLUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 54,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 636,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                        {
                                            "id": "NP_PLUS",
                                            "type": "CDS",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 562,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                    ],
                                }
                            ],
                        }
                    ],
                },
                "sequence": {
                    "seq": "auggggucacagugguuauaaagaguuauacccugaagaauuugaaacagacaguagugaucagcaagauauuaccaacgggaagaaaacaucuccccagguaaagucaucuacccaugaaucccgcaaacacaagaagucaaagaaaucccacaaaaaaaagcagaaaaaaaggucacacaaaaaacagaagaaaagcaaaaaggaagccacagauauaacagcagauuccucgagugaguucucagaagaaacuggggcuucugguacaaggaaagggaaacaaccacauaaacgcaagaaaaaauccaggaaaaagucucucaaaaaaccugcuuuauucuuagaggcagaaaguaacacuucacauucagaugauucagcauccagcaguucugaggaaagugaggaaagagacacuaagaaaaccaaaaggaaaaagagagagaaaaaagcccauaccucuguagccaacaaugaaauacaggagaggacaaacaaacgcacaaauuggaaaguagcuacagaugaaaggucugcugagagcucagaggaugacuaaaugggaaacacuuuuguuuuccacaugacuguggauauuuacaguucuuacuccuugugguuuugccagugacu"
                },
            },
        ),
        (
            "TEST_REF",
            "NR_PLUS",
            {
                "annotations": {
                    "qualifiers": {
                        "name": "TEST_REF",
                        "chromosome": "TEST_CHR",
                        "mol_type": "genomic DNA",
                    },
                    "id": "TEST_REF",
                    "type": "record",
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 0},
                        "end": {"type": "point", "position": 503},
                    },
                    "features": [
                        {
                            "qualifiers": {"name": "PLUS"},
                            "id": "PLUS",
                            "type": "gene",
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": 0},
                                "end": {"type": "point", "position": 503},
                                "strand": 1,
                            },
                            "features": [
                                {
                                    "id": "NR_PLUS",
                                    "type": "ncRNA",
                                    "location": {
                                        "type": "range",
                                        "start": {"type": "point", "position": 0},
                                        "end": {"type": "point", "position": 503},
                                        "strand": 1,
                                    },
                                    "features": [
                                        {
                                            "id": "exon-NR_PLUS-1",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 101,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                        {
                                            "id": "exon-NR_PLUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 101,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 202,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                        {
                                            "id": "exon-NR_PLUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 202,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 503,
                                                },
                                                "strand": 1,
                                            },
                                        },
                                    ],
                                }
                            ],
                        }
                    ],
                },
                "sequence": {
                    "seq": "uugcuuuaacaauacaugugaugugucauauuacagauggggucacagugguuauaaagaguuauacccugaagaauuugaaacagacagguaaggaaaaugaaacaauaggugauuuaguuaggaggugaucccuguuuucccauucucugguugauguuuggcaugucuguaagcauuuugguuuuuauauauaguauucuuuuaacuuuucuuauuaguagugaucagcaagauauuaccaacgggaagaaaacaucuccccagguaaagucaucuacccaugaaucccgcaaacacaagaagucaaagaaaucccacaaaaaaaagcagaaaaaaaggucacacaaaaaacagaagaaaagcaaaaaggaagccacagauauaacagcagauuccucgagugaguucucagaagaaacuggggcuucugguacaaggaaagggaaacaaccacauaaacgcaagaaaaaauccaggaaaaagucucucaaaaaaccugc"
                },
            },
        ),
        (
            "TEST_REF",
            "NM_MINUS",
            {
                "annotations": {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 0},
                        "end": {"type": "point", "position": 954},
                    },
                    "qualifiers": {
                        "name": "TEST_REF",
                        "chromosome": "TEST_CHR",
                        "mol_type": "genomic DNA",
                    },
                    "type": "record",
                    "id": "TEST_REF",
                    "features": [
                        {
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": 0},
                                "end": {"type": "point", "position": 954},
                                "strand": -1,
                            },
                            "qualifiers": {"name": "MINUS"},
                            "type": "gene",
                            "id": "MINUS",
                            "features": [
                                {
                                    "id": "NM_MINUS",
                                    "type": "mRNA",
                                    "location": {
                                        "type": "range",
                                        "start": {"type": "point", "position": 0},
                                        "end": {"type": "point", "position": 954},
                                        "strand": -1,
                                    },
                                    "features": [
                                        {
                                            "id": "exon-NM_MINUS-1",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 603,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 954,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NM_MINUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 502,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 603,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NM_MINUS-3",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 201,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 502,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NM_MINUS-4",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 201,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "NP_MINUS",
                                            "type": "CDS",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 125,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 854,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                    ],
                                }
                            ],
                        }
                    ],
                },
                "sequence": {
                    "seq": "uguacuccaaagucuuucuaauguugcuuuaauuuccaaaaauguaugcauugcuuuaacaauacaugugaugugucauauuacagauggggucacagugguuauaaagaguuauacccugaagaauuugaaacagacagguaaggaaaauaggcuuacugaaagaaacuaagaugguacaaaaucuguauuauaaauugagguugauguuuggcaugucuguaagcauuuugguuuuuauauauaguauuccauacaguaacucaguauggcagcuuagaauuuuuaccuucauuuuaaagaugaggaaacaaaaacucaaugagaauauuaaaguguuaaaguauacauuaaagugcuuauuuaaaauucagauguuaaccucaauuuuuuaaucuagaaugcaaaauauuaaaauaauacgcuuuuuuuuuacauaaaagcuucuauuuuuuaacuuuucuuauuaguagugaucagcaagauauuaccaacgggaagaagugaguucucagaagaaacuggggcuucugguacaaggaaagggaaacaaccacauaaacgcaagaaaaaauccaggaaaaagucucucaaaaaaccugcgaaaaagagagagaaaaaagcccauaccucuguagccaacaaugaaauacaggagaggacaaacaaacgcacaaauuggaaaguagcuacagaugaaaggucugcugagagcucagaggaugacuaaaugggaaacacuuuuguuuuccacaugacuguggauauuuacaguucuuacuccuugugguuuugccagugacucuuguucagcacggggccugaggucagagcugucuugugccaucugucauuucugacagacgucuugucuucuauuuuggcguuaagcuugauccccuuuucuuguuaaaagggaaucugguauuuuguuaugaagguuucuugaagaga"
                },
            },
        ),
        (
            "TEST_REF",
            "NR_MINUS",
            {
                "annotations": {
                    "location": {
                        "type": "range",
                        "start": {"type": "point", "position": 0},
                        "end": {"type": "point", "position": 993},
                    },
                    "id": "TEST_REF",
                    "qualifiers": {
                        "name": "TEST_REF",
                        "chromosome": "TEST_CHR",
                        "mol_type": "genomic DNA",
                    },
                    "type": "record",
                    "features": [
                        {
                            "location": {
                                "type": "range",
                                "start": {"type": "point", "position": 0},
                                "end": {"type": "point", "position": 993},
                                "strand": -1,
                            },
                            "id": "MINUS",
                            "qualifiers": {"name": "MINUS"},
                            "type": "gene",
                            "features": [
                                {
                                    "id": "NR_MINUS",
                                    "type": "ncRNA",
                                    "location": {
                                        "type": "range",
                                        "start": {"type": "point", "position": 0},
                                        "end": {"type": "point", "position": 993},
                                        "strand": -1,
                                    },
                                    "features": [
                                        {
                                            "id": "exon-NR_MINUS-1",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 542,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 993,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NR_MINUS-2",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 191,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 542,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                        {
                                            "id": "exon-NR_MINUS-3",
                                            "type": "exon",
                                            "location": {
                                                "type": "range",
                                                "start": {
                                                    "type": "point",
                                                    "position": 0,
                                                },
                                                "end": {
                                                    "type": "point",
                                                    "position": 191,
                                                },
                                                "strand": -1,
                                            },
                                        },
                                    ],
                                }
                            ],
                        }
                    ],
                },
                "sequence": {
                    "seq": "agucuuucuaauguugcuuuaauuuccaaaaauguaugcauugcuuuaacaauacaugugaugugucauauuacagauggggucacagugguuauaaagaguuauacccugaagaauuugaaacagacagguaaggaaaauaggcuuacugaaagaaacuaagaugguacaaaaucuguauuauaaauugagguugauguuuggcaugucuguaagcauuuugguuuuuauauauaguauuccauacaguaacucaguauggcagcuuagaauuuuuaccuucauuuuaaagaugaggaaacaaaaacucaaugagaauauuaaaguguuaaaguauacauuaaagugcuuauuuaaaauucagauguuaaccucaauuuuuuaaucuagaaugcaaaauauuaaaauaauacgcuuuuuuuuuacauaaaagcuucuauuuuuuaacuuuucuuauuaguagugaucagcaagauauuaccaacgggaagaaaacaucuccccagguaaagucaucuacccaugaaucccgcaaacacaagcuuuauucuuagaggcagaaaguaacacuucacauucagaugauucagcauccagcaguucugaggaaagugaggaaagagacacuaagaaaaccaaaaggaaaaagagagagaaaaaagcccauaccucuguagccaacaaugaaauacaggagaggacaaacaaacgcacaaauuggaaaguagcuacagaugaaaggucugcugagagcucagaggaugacuaaaugggaaacacuuuuguuuuccacaugacuguggauauuuacaguucuuacuccuugugguuuugccagugacucuuguucagcacggggccugaggucagagcugucuugugccaucugucauuucugacagacgucuugucuucuauuuuggcguuaagcuugauccccuuuucuuguuaaaagggaaucugguauuuuguuaugaagguuucuugaagaga"
                },
            },
        ),
    ],
)
def test_to_rna_reference_model(r_id, s_id, expected):
    model = retrieve_reference(r_id)[0]
    rna_model = to_rna_reference_model(model, s_id, transcribe=True)
    assert rna_model == expected


def _location(s, e, strand=None):
    location = {
        "type": "range",
        "start": {"type": "point", "position": s},
        "end": {"type": "point", "position": e},
    }
    if strand:
        location["strand"] = strand
    return location


def _variant(s, e, ins=None):
    v = {
        "location": _location(s, e),
        "type": "deletion_insertion",
        "source": "reference",
    }
    if ins:
        v["inserted"] = [{"sequence": ins, "source": "description"}]
    return v


TESTS_VARIANTS = [
    # input, trimmed, rna
    (  # same intron
        [_variant(10, 20, "T")],
        [],
        [],
    ),
    (  # intron exon with insertion
        [_variant(120, 141, "A")],
        [],
        [],
    ),
    (  # same exon
        [_variant(150, 180, "A")],
        [_variant(150, 180, "A")],
        [_variant(15, 45, "a")],
    ),
    (  # intron intron with insertion
        [_variant(10, 200, "A")],
        [],
        [],
    ),
    (  # intron intron without insertion
        [_variant(10, 200)],
        [_variant(135, 189)],
        [_variant(0, 54)],
    ),
    (
        [_variant(10, 189, "A")],
        [],
        [],
    ),
    (  # intron intron over exon with insertion
        [_variant(10, 1200, "A")],
        [],
        [],
    ),
    (  # intron intron over exon without insertion
        [_variant(10, 1210)],
        [_variant(135, 1200)],
        [_variant(0, 636)],
    ),
    (  # exon - exon over intron without insertion
        [_variant(150, 1100)],
        [_variant(150, 1100)],
        [_variant(15, 536)],
    ),
]


@pytest.mark.parametrize(
    "variants, expected",
    [(t_v[0], t_v[1]) for t_v in TESTS_VARIANTS],
)
def test_trim_to_exons(variants, expected):
    sequences = {
        "TEST_REF": retrieve_reference("TEST_REF")[0],
        "reference": retrieve_reference("TEST_REF")[0],
    }
    selector_model = get_internal_selector_model(
        retrieve_reference("TEST_REF")[0]["annotations"], "NM_PLUS"
    )

    # exons: [135, 189, 618, 1200]
    assert _trim_to_exons(variants, selector_model["exon"], sequences) == expected


@pytest.mark.parametrize(
    "variants, expected",
    [(t_v[0], t_v[2]) for t_v in TESTS_VARIANTS],
)
def test_to_rna_variants(variants, expected):
    sequences = {
        "TEST_REF": retrieve_reference("TEST_REF")[0],
        "reference": retrieve_reference("TEST_REF")[0],
    }
    selector_model = get_internal_selector_model(
        retrieve_reference("TEST_REF")[0]["annotations"], "NM_PLUS"
    )
    # exons: [135, 189, 618, 1200]
    assert to_rna_variants(variants, sequences, selector_model) == expected


@pytest.mark.parametrize(
    "location, exons, location_type",
    [
        ("135_135", [(135, 189), (618, 1200)], "same exon"),
        ("135_136", [(135, 189), (618, 1200)], "same exon"),
        ("135_188", [(135, 189), (618, 1200)], "same exon"),
        ("135_189", [(135, 189), (618, 1200)], "same exon"),
        ("135_619", [(135, 189), (618, 1200)], "exon exon"),
        ("136_188", [(135, 189), (618, 1200)], "same exon"),
        ("136_189", [(135, 189), (618, 1200)], "same exon"),
        ("136_619", [(135, 189), (618, 1200)], "exon exon"),
        ("188_188", [(135, 189), (618, 1200)], "same exon"),
        ("188_189", [(135, 189), (618, 1200)], "same exon"),
        ("618_618", [(135, 189), (618, 1200)], "same exon"),
        ("618_619", [(135, 189), (618, 1200)], "same exon"),
        ("618_1199", [(135, 189), (618, 1200)], "same exon"),
        ("618_1200", [(135, 189), (618, 1200)], "same exon"),
        ("1199_1199", [(135, 189), (618, 1200)], "same exon"),
        ("1199_1200", [(135, 189), (618, 1200)], "same exon"),
        ("188_188", [(135, 189), (189, 1200)], "same exon"),
        ("188_189", [(135, 189), (189, 1200)], "same exon"),
        ("189_189", [(135, 189), (189, 1200)], "same exon"),
        ("188_190", [(135, 189), (189, 1200)], "exon exon"),
    ],
)
def test_get_location_type(location, exons, location_type):
    location = to_model(location, "location")
    assert get_location_type(location, exons) == location_type


@pytest.mark.parametrize(
    "position, exons, position_type",
    [
        (1, [(10, 20), (30, 40), (40, 50)], (0, 0)),
        (2, [(10, 20), (30, 40), (40, 50)], (0, 0)),
        (3, [(10, 20), (30, 40), (40, 50)], (0, -2)),
        (4, [(10, 20), (30, 40), (40, 50)], (0, -2)),
        (5, [(10, 20), (30, 40), (40, 50)], (0, -2)),
        (6, [(10, 20), (30, 40), (40, 50)], (0, -2)),
        (7, [(10, 20), (30, 40), (40, 50)], (0, -2)),
        (8, [(10, 20), (30, 40), (40, 50)], (0, -1)),
        (9, [(10, 20), (30, 40), (40, 50)], (0, -1)),
        (10, [(10, 20), (30, 40), (40, 50)], (1, 0)),
        (11, [(10, 20), (30, 40), (40, 50)], (1, 0)),
        (12, [(10, 20), (30, 40), (40, 50)], (1, 0)),
        (19, [(10, 20), (30, 40), (40, 50)], (1, 0)),
        (20, [(10, 20), (30, 40), (40, 50)], (2, 1)),
        (21, [(10, 20), (30, 40), (40, 50)], (2, 1)),
        (22, [(10, 20), (30, 40), (40, 50)], (2, 2)),
        (23, [(10, 20), (30, 40), (40, 50)], (2, 2)),
        (24, [(10, 20), (30, 40), (40, 50)], (2, 2)),
        (25, [(10, 20), (30, 40), (40, 50)], (2, -2)),
        (26, [(10, 20), (30, 40), (40, 50)], (2, -2)),
        (27, [(10, 20), (30, 40), (40, 50)], (2, -2)),
        (28, [(10, 20), (30, 40), (40, 50)], (2, -1)),
        (29, [(10, 20), (30, 40), (40, 50)], (2, -1)),
        (30, [(10, 20), (30, 40), (40, 50)], (3, 0)),
        (31, [(10, 20), (30, 40), (40, 50)], (3, 0)),
        (39, [(10, 20), (30, 40), (40, 50)], (3, 0)),
        (40, [(10, 20), (30, 40), (40, 50)], (5, 0)),
        (41, [(10, 20), (30, 40), (40, 50)], (5, 0)),
        (49, [(10, 20), (30, 40), (40, 50)], (5, 0)),
        (50, [(10, 20), (30, 40), (40, 50)], (6, 1)),
        (51, [(10, 20), (30, 40), (40, 50)], (6, 1)),
        (52, [(10, 20), (30, 40), (40, 50)], (6, 2)),
        (57, [(10, 20), (30, 40), (40, 50)], (6, 0)),
        (60, [(10, 20), (30, 40), (40, 50)], (6, 0)),
        # ---
        (0, [(0, 10), (10, 20), (20, 30)], (1, 0)),
        (9, [(0, 10), (10, 20), (20, 30)], (1, 0)),
        (10, [(0, 10), (10, 20), (20, 30)], (3, 0)),
    ],
)
def test_get_position_type(position, exons, position_type):
    assert get_position_type(position, exons, 2, 5) == position_type


@pytest.mark.parametrize(
    "variant, error, info", [(_variant(10, 200), [], ["ISPLICESITEREMOVED"])]
)
def test_splice_sites(variant, error, info):
    sequences = {
        "TEST_REF": retrieve_reference("TEST_REF")[0],
        "reference": retrieve_reference("TEST_REF")[0],
    }
    selector_model = get_internal_selector_model(
        retrieve_reference("TEST_REF")[0]["annotations"], "NM_PLUS"
    )
    # exons: [135, 189, 618, 1200]
    errors, infos = splice_sites([variant], sequences, selector_model)

    if error:
        assert errors[0]["code"] == error[0]
    if info:
        assert infos[0]["code"] == info[0]
