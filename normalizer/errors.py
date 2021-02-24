from .description_model import (
    location_to_description,
    point_to_description,
    variant_to_description,
)


def mismatch(input_description, model_description):
    return {
        "code": "EMISMATCH",
        "details": "Model description {} different than the input description {}.".format(
            model_description, input_description
        ),
    }


def reference_not_retrieved(reference_id, path):
    return {
        "code": "ERETR",
        "details": "Reference {} could not be retrieved.".format(reference_id),
        "paths": [path],
    }


def no_selector_found(reference_id, selector_id, path):
    return {
        "code": "ENOSELECTORFOUND",
        "details": "No {} selector found in reference {}.".format(
            selector_id, reference_id
        ),
        "paths": [path],
    }


def selector_options(selector_id, selector_type, options, path):
    return {
        "code": "ESELECTOROPTIONS",
        "details": "{} selector identified as {}.".format(selector_id, selector_type),
        "options": options,
        "paths": [path],
    }


def no_coordinate_system(path):
    return {
        "code": "ENOCOORDINATESYSTEM",
        "details": "A coordinate system is required.",
        "paths": [path],
    }


def coordinate_system_mismatch(
    coordinate_system, mismatch_id, mismatch_coordinate_system, path
):
    return {
        "code": "ECOORDINATESYSTEMMISMATCH",
        "details": "Coordinate system {} does not match with {} "
        "{} coordinate system. ".format(
            coordinate_system, mismatch_id, mismatch_coordinate_system
        ),
        "paths": [path],
    }


def offset(location, path):
    return {
        "code": "EOFFSET",
        "details": "Offsets, as in `{}', are not allowed with the 'g' coordinate system.".format(
            location_to_description(location)
        ),
        "paths": [path],
    }


def outside_cds(location, path):
    if location["outside_cds"] == "upstream":
        d_in = "-"
    elif location["outside_cds"] == "downstream":
        d_in = "*"
    else:
        d_in = ""
    return {
        "code": "EOUTSIDECDS",
        "details": "Outside CDS specifics, as '{}' in '{}', are not allowed with the 'g' coordinate system.".format(
            d_in, location_to_description(location)
        ),
        "paths": [path],
    }


def out_of_boundary_lesser(position, path):
    return {
        "code": "EOUTOFBOUNDARY",
        "details": "Position {} is lesser than 1.".format(
            point_to_description(position)
        ),
        "paths": [path],
    }


def out_of_boundary_greater(point, sequence_length, path):
    return {
        "code": "EOUTOFBOUNDARY",
        "details": "Position {} is greater than the sequence {} length.".format(
            point_to_description(point), sequence_length
        ),
        "paths": [path],
    }


def range_reversed(location, path):
    return {
        "code": "ERANGEREVERSED",
        "details": "Range start position greater than the end position.".format(
            location_to_description(location)
        ),
        "paths": [path],
    }


def insertion_range(location, path):
    return {
        "code": "EINSERTIONRANGE",
        "details": "Range positions {} not consecutive in insertion location.".format(
            location_to_description(location)
        ),
        "paths": [path],
    }


def insertion_location_not_range(point, path):
    return {
        "code": "EINSERTIONLOCATION",
        "details": "Insertion location {} is not range.".format(
            point_to_description(point)
        ),
        "paths": [path],
    }


def repeat_reference_sequence_length(path):
    return {
        "code": "EREPEATREFERENCELENGTH",
        "details": "Reference sequence length not a multiple of the inserted sequence length.",
        "paths": [path],
    }


def repeat_sequences_mismatch(reference_sequence, repeat_sequence, path):
    return {
        "code": "EREPEATMISMATCH",
        "details": "Reference sequence {} does not contain the repeat sequence {}.".format(
            reference_sequence, repeat_sequence
        ),
        "paths": [path],
    }


def length_mismatch(reference_length, deleted_length, path):
    return {
        "code": "ELENGTHMISMATCH",
        "details": "The length {} differs from that of the range {}.".format(
            deleted_length, reference_length
        ),
        "paths": [path],
    }


def sequence_mismatch(reference_sequence, deleted_sequence, path):
    return {
        "code": "ESEQUENCEMISMATCH",
        "details": "{} not found in the reference sequence, found {} instead.".format(
            deleted_sequence, reference_sequence
        ),
        "paths": [path],
    }


def no_dna(sequence, path):
    return {
        "code": "ENODNA",
        "details": "Sequence '{}' is not a DNA sequence.".format(sequence),
        "paths": [path],
    }


def repeat_not_supported(variant, path):
    return {
        "code": "EREPEATUNSUPPORTED",
        "details": "Repeat variant {} not supported.".format(
            variant_to_description(variant)
        ),
        "paths": [path],
    }


def overlap():
    return {"code": "EOVERLAP", "details": "Variant locations overlap."}


def syntax_uc(e):
    return dict(
        {"code": "ESYNTAXUC", "details": "Unexpected character."}, **e.serialize()
    )


def syntax_ueof(e):
    return dict(
        {"code": "ESYNTAXUEOF", "details": "Unexpected end of input."}, **e.serialize()
    )


def position_syntax(details, e):
    return dict({"code": "EPOSITIONSYNTAX", "details": details}, **e.serialize())


def position_invalid():
    return {"code": "EPOSITIONINVALID", "details": "Position must be string"}


def no_inputs():
    return {"code": "ENOINPUTS"}


def no_inputs_other():
    return {"code": "ENOINPUTSOTHER"}


def no_to_selector():
    return {"code": "ENOTOSELECTOR"}
