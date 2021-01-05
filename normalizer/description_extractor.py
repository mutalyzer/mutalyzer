import extractor

from .converter import de_to_hgvs
from .converter.to_hgvs_coordinates import (
    crossmap_to_hgvs_setup,
    locations_to_hgvs_locations,
)
from .converter.to_hgvs_indexing import variants_to_internal_indexing
from .description_model import variants_to_description


def description_extractor(reference, observed):
    de_variants = extractor.describe_dna(reference, observed)
    de_hgvs_variants = de_to_hgvs(
        de_variants, {"reference": reference, "observed": observed}
    )
    hgvs_indexing_variants = variants_to_internal_indexing(de_hgvs_variants)
    crossmap = crossmap_to_hgvs_setup("g")
    hgvs_variants = locations_to_hgvs_locations({"variants": hgvs_indexing_variants}, crossmap)
    normalized_description = variants_to_description(hgvs_variants["variants"])
    return normalized_description
