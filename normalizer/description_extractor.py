import extractor

from .converter import de_to_hgvs, to_hgvs
from .description import variants_to_description


def description_extractor(reference, observed):
    de_variants = extractor.describe_dna(reference, observed)
    de_hgvs_variants = de_to_hgvs(
        de_variants, {"reference": reference, "observed": observed}
    )
    crossmap = to_hgvs.crossmap_coordinate_to_genomic_setup()
    hgvs_variants = to_hgvs.variants_to_hgvs(de_hgvs_variants, crossmap)
    normalized_description = variants_to_description(hgvs_variants)
    return normalized_description
