# from algebra.extractor import extract_sequence, to_hgvs
from algebra import LCSgraph, Variant
from .algebra import to_hgvs

def description_extractor(reference, observed):
    graph = LCSgraph.from_variants(reference, [Variant(0, len(reference), observed)])
    variants = graph.canonical()
    # variants, _ = extract_sequence(reference, observed)
    return "[" + ";".join([to_hgvs(v) for v in variants]) + "]"
