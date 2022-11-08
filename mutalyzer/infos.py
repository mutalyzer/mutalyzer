from .description_model import point_to_description, variant_to_description


def corrected_reference_id(original_id, corrected_id, path):
    return {
        "code": "ICORRECTEDREFERENCEID",
        "details": "Reference {} was retrieved instead of {}.".format(
            corrected_id, original_id
        ),
        "paths": [path],
    }


def corrected_lrg_reference(original_id, lrg_model, path):
    return {
        "code": "ICORRECTEDLRGREFERENCE",
        "details": "Reference {} was identified as LRG with ID {} and selector {}.".format(
            original_id, lrg_model["id"], lrg_model["selector"]["id"]
        ),
        "paths": [path],
    }


def corrected_selector_id(original_id, corrected_id, correction_source, path):
    return {
        "code": "ICORRECTEDSELECTORID",
        "details": "Selector {} was corrected to {} from {}.".format(
            original_id, corrected_id, correction_source
        ),
        "paths": [path],
    }


def corrected_coordinate_system(coordinate_system, correction_source, path):
    return {
        "code": "ICORRECTEDCOORDINATESYSTEM",
        "details": "Coordinate system corrected to {} from {}.".format(
            coordinate_system, correction_source
        ),
        "paths": [path],
    }


def to_coordinate_system_from_reference(coordinate_system):
    return {
        "code": "ITOSELECTORFROMREFERENCE",
        "details": "To coordinate system identified as {} "
        "from the selector molecule type.".format(coordinate_system),
    }


def corrected_variant_type(original_type, corrected_type):
    return {
        "code": "ICORRECTEDVARIANTTYPE",
        "details": "Variant corrected from '{}' to '{}'".format(
            original_type, corrected_type
        ),
    }


def corrected_point(original, corrected, path):
    return {
        "code": "ICORRECTEDPOINT",
        "details": "Point corrected from '{}' to '{}'".format(
            point_to_description(original), point_to_description(corrected)
        ),
        "paths": [path],
    }


def from_to_selector_equal():
    return {
        "code": "IFROMTOSELECTORSEQUAL",
        "details": "From and to coordinate systems are equal.",
    }


def sorted_variants():
    return {
        "code": "ISORTEDVARIANTS",
        "details": "Variants were sorted according to their locations.",
    }


def corrected_sequence(original, corrected):
    return {
        "code": "CORRECTEDSEQUENCE",
        "details": f'Sequence "{original}" corrected to "{corrected}".',
    }


def variant_discarded(path):
    return {
        "code": "IVARIANTDISCARDED",
        "details": "Variant discarded.",
        "paths": [path],
    }


def splice_site_removed(path):
    return {
        "code": "ISPLICESITEREMOVED",
        "details": "Splice(s) sites removed.",
        "paths": [path],
    }


def whole_transcript_exon(reference_id, selector_id, path):
    return {
        "code": "IWHOLETRANSCRIPTEXON",
        "details": f"No exon features found in the '{reference_id}' reference"
        f" sequence for '{selector_id}'. The entire transcript was"
        f" assumed as one exon.",
        "paths": [path],
    }


def insertions_same_location(variants, paths):
    plural = "s" if len(variants) > 1 else ""
    variants = ", ".join([variant_to_description(v) for v in variants])
    return {
        "code": "IINSERTIONSSAMELOCATION",
        "details": f"The following insertion{plural} {variants} are at the same location.",
        "paths": paths,
    }


def mrna_genomic_tip():
    return {
        "code": "IMRNAGENOMICTIP",
        "details": "An 'mRNA' sequence was used with the 'c.' coordinate system."
        " Make use of a genomic reference sequence if the experiment "
        "performed involved measured DNA.",
    }


def mrna_genomic_difference(mrna_id, genomic_id):
    return {
        "code": "IMRNAGENOMICDIFFERENCE",
        "details": f"There are differences between the mRNA sequence of {mrna_id} and the genomic sequence of {genomic_id}.",
    }


def no_selector(reference_id, selector_id):
    return {
        "code": "INOSELECTOR",
        "details": "No {} selector found in {}.".format(selector_id, reference_id),
    }


def other_versions(reference_id, selector_id, other_versions):
    multiple = "s" if len(other_versions) > 1 else ""
    return {
        "code": "IOTHERVERSIONS",
        "details": f"Selector id{multiple} {', '.join(other_versions)} found in {reference_id} instead of {selector_id}.",
    }


def assembly_chromosome_to_id(assembly_id, chromosome_number, chromosome_id):
    return {
        "code": "ASSEMBLY_CHROMOSOME_TO_ID",
        "details": f"Assembly {assembly_id} and {chromosome_number} corrected to {chromosome_id}.",
    }
