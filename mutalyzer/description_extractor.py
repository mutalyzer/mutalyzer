from algebra.extractor import extract_sequence, to_hgvs


def description_extractor(reference, observed):
    variants, _ = extract_sequence(reference, observed)
    return to_hgvs(variants, reference)
