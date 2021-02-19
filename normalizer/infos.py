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
