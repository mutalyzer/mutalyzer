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
