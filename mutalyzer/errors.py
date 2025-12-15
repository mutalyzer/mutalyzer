from .description_model import (
    location_to_description,
    point_to_description,
    variant_to_description,
    point_outside_cds,
    point_position,
)


def mismatch(input_description, model_description):
    return {
        "code": "EMISMATCH",
        "details": f"Model description {model_description} differs from the input description {input_description}.",
    }

def gene_as_reference_id(gene_name, chr_ids, options, path):
    # SDHD:c.274G>T
    return {
        "code": "EGENEASREFERENCEID",
        "details": f"`{gene_name}` is an invalid reference identifier; it appears to be a gene name.",
        "gene": gene_name,
        "chr_ids": chr_ids,
        "options": options,
        "paths": [path],
    }

def selector_options(selector_id, selector_type, options, path):
    # Gene name only: NG_007485.1(CDKN2A):n.204_205insATC
    # Gene name with some old format: NG_007485.1(CDKN2A_v001):n.204_205insATC
    # HGNC id: NG_007485.1(1787):n.204_205insATC
    # Check if it should be merged with EGENEASREFERENCEID.
    # For a transcript accession without a version it could be raised as well.
    return {
        "code": "ESELECTOROPTIONS",
        "details": f"`{selector_id}` is an invalid selector; it was identified as a `{selector_type}` selector.",
        "options": options,
        "paths": [path],
    }

def reference_not_retrieved(reference_id, path):
    # NO_REF:g.100del
    return {
        "code": "ERETR",
        "details": f"`{reference_id}` could not be retrieved; it may be an invalid reference identifier.",
        "paths": [path],
    }


def no_selector_found(reference_id, selector_id, path):
    # NG_012337.1(DUMMYACCNO_9999.9):c.12_13insGATC
    # NG_012337.3(NM_003002.1):c.274G>T -> we may suggest that a newer version (NM_003002.4) is available.
    return {
        "code": "ENOSELECTORFOUND",
        "details": f"`{selector_id}` is not annotated in `{reference_id}`.",
        "paths": [path],
    }


def coordinate_system_invalid(coordinate_system, path):
    # NG_012337.3:a.274G>T
    return {
        "code": "ECOORDINATESYSTEMINVALID",
        "details": f"`{coordinate_system}.` is an invalid coordinate system.",
        "paths": [path],
    }


def no_coordinate_system(path):
    # It was not possible to identify the coordinate system.
    return {
        "code": "ENOCOORDINATESYSTEM",
        "details": "A coordinate system (e.g., `g.`, `c.`, `r.`, ...) is required.",
        "paths": [path],
    }


def coordinate_system_mismatch(coordinate_system, mismatch_id, mismatch_coordinate_system, path):
    # NR_038420.1:c.10del
    # Check if it can be merged with ECOORDINATESYSTEMINVALID
    # Check if there could be multiple options (NM with c. or r.).
    return {
        "code": "ECOORDINATESYSTEMMISMATCH",
        "details": f"`{coordinate_system}.` is an invalid coordinate system in the context of {mismatch_id}; "
                   f"`{mismatch_coordinate_system}.` is expected, for example.",
        "paths": [path],
    }


def offset(location, path):
    # NG_012337.1:g.7125+1G>T
    # Check the ` consistency.
    return {
        "code": "EOFFSET",
        "details": f"`{location_to_description(location)}` contains an offset, which is not allowed with the `g.` "
                   f"coordinate system.",
        "paths": [path],
    }


def offset_direction(point, path):
    # NG_012337.3(NM_003002.4):c.52-20del
    sign = "-" if point["offset"]["value"] < 0 else "+"
    return {
        "code": "EOFFSETDIRECTION",
        "details": f"`{location_to_description(point)}` is invalid; the sign `{sign}` of the offset may be wrong "
                   f"or the wrong exon boundary was used.",
        "paths": [path],
    }


def exon_boundary(point, path):
    return {
        "code": "EEXONBOUNDARY",
        "details": f"`{location_to_description(point)}` is invalid; `{point_position(point)}` is not an exon boundary.",
        "paths": [path],
    }


def outside_cds(point, path):
    # NG_007485.1:g.-1del
    return {
        "code": "EOUTSIDECDS",
        "details": f"`{location_to_description(point)}` is invalid; "
                   f"outside CDS specifics such as `{point_outside_cds(point)}` are not allowed with the `g.` "
                   f"coordinate system.",
        "paths": [path],
    }


def intronic(positions, options=None, assemblies=None):
    # NM_003002.4:c.50+10del
    details_parts = []
    for ref_id, location_descs in positions.items():
        if len(location_descs) > 1:
            locations = ", ".join(f"`{d}`" for d in location_descs)
            details_parts.append(
                f"Intronic positions {locations} were used with a non intronic reference sequence `{ref_id}`."
            )
        else:
            details_parts.append(
                f"Intronic position `{location_descs[0]}` was used with a non intronic reference sequence `{ref_id}`."
            )

    output = {
        "code": "EINTRONIC",
        "positions": positions,
        "details": " ".join(details_parts)
    }

    if options:
        output["options"] = options
    if assemblies:
        output["assemblies"] = assemblies

    return output

def intronic_rna(location, path):
    # "NG_012337.3(NM_003002.4):r.274+10G>T
    return {
        "code": "EINTRONICRNA",
        "details": f"Intronic position `{location_to_description(location)}` was given for an RNA description.",
        "paths": [path],
    }


def out_of_boundary_lesser(position, shift, path):
    # NG_007485.1(NM_000077.4):c.0del
    plural = "s" if shift > 1 else ""
    return {
        "code": "EOUTOFBOUNDARY",
        "details": f"Position `{point_to_description(position)}` is invalid; it is {shift} nucleotide{plural} before "
                   f"the start of the sequence.",
        "paths": [path],
    }


def out_of_boundary_greater(point, shift, sequence_length, path):
    # NG_007485.1(NM_000077.4):c.40000del
    plural = "s" if shift > 1 else ""
    return {
        "code": "EOUTOFBOUNDARY",
        "details": f"Position `{point_to_description(point)}` is invalid; it is {shift} nucleotide{plural} after the "
                   f"end of the sequence (sequence length is {sequence_length} nucleotides).",
        "paths": [path],
    }


def range_reversed(location, path):
    # "NG_007485.1:g.4_3del
    return {
        "code": "ERANGEREVERSED",
        "details": f"Variant range is invalid; start position is greater than the end position in "
                   f"`{location_to_description(location)}`.",
        "paths": [path],
    }


def insertion_range(location, path):
    # NG_007485.1:g.40_42insT
    return {
        "code": "EINSERTIONRANGE",
        "details": f"Insertion variant is invalid; range positions `{location_to_description(location)}` are not "
                   f"consecutive.",
        "paths": [path],
    }


def repeat_reference_sequence_length(len_segment, len_unit, path):
    # NG_012337.1(NM_003002.2):c.100_102AA[4]
    return {
        "code": "EREPEATREFERENCELENGTH",
        "details": f"Repeat variant is invalid; reference segment length `{len_segment}` is not a multiple of the "
                   f"repeat unit length `{len_unit}`.",
        "paths": [path],
    }


def repeat_sequences_mismatch(reference_sequence, repeat_sequence, path):
    # NG_012337.1(NM_003002.2):c.100_102AAT[4]
    return {
        "code": "EREPEATMISMATCH",
        "details": f"Repeat variant is invalid; reference sequence `{reference_sequence}` does not contain the "
                   f"repeat sequence `{repeat_sequence}`.",
        "paths": [path],
    }


def length_mismatch(reference_length, deleted_length, path):
    # NG_012337.1(NM_003002.2):c.274del3
    return {
        "code": "ELENGTHMISMATCH",
        "details": f"Deleted length `{deleted_length}` differs from the reference range length `{reference_length}`.",
        "paths": [path],
    }


def sequence_mismatch(reference_sequence, deleted_sequence, path):
    # NG_012337.1(NM_003002.2):c.274A>T
    return {
        "code": "ESEQUENCEMISMATCH",
        "details": f"`{deleted_sequence}` was not found in the reference sequence; `{reference_sequence}` was found "
                   f"instead.",
        "paths": [path],
    }


def amino_acid_mismatch(description_aa, reference_aa, path):
    # NM_005410.4:p.D59delinsA
    return {
        "code": "EAMINOACIDMISMATCH",
        "details": f"Amino acid sequence `{description_aa}` was not found in the reference sequence; `{reference_aa}` "
                   f"was found instead.",
        "paths": [path],
    }


def no_dna(sequence, path):
    # NM_000143.3:c.45delU"
    return {
        "code": "ENODNA",
        "details": f"Sequence `{sequence}` is invalid; it is not a DNA sequence.",
        "paths": [path],
    }


def no_rna(sequence, path):
    # NM_003002.2:r.277t>u
    return {
        "code": "ENORNA",
        "details": f"Sequence `{sequence}` is invalid; it is not an RNA sequence.",
        "paths": [path],
    }


def repeat_not_supported(variant, path):
    # LRG_1t1:52_153CAG[21]CAA[1]CAG[1]CCG[1]CCA[1]CCG[7]CCT[2]
    return {
        "code": "EREPEATUNSUPPORTED",
        "details": f"Repeat variant `{variant_to_description(variant)}` is not supported.",
        "paths": [path],
    }


def variant_not_supported(variant, variant_type, path):
    # NG_007485.1(NP_000068.1):p.(Met54Ilefs*66)
    return {
        "code": "EVARIANTNOTSUPPORTED",
        "details": f"Variant `{variant_to_description(variant)}` of type `{variant_type}` is not supported.",
        "paths": [path],
    }


def uncertain():
    # Improve support for some uncertain variants.
    return {"code": "EUNCERTAIN", "details": "Cannot proceed; uncertainties are present in one or more locations."}


def inserted_length():
    # NG_012337.1:g.20_21ins[123;30_40]
    # Check if we should still support lengths as such.
    return {
        "code": "EINSERTEDLENGTH",
        "details": "Inserted part length is not supported.",
    }


def overlap():
    return {"code": "EOVERLAP", "details": "Variants overlap; overlapping variant locations are not allowed."}


def syntax_uc(e):
    return dict({"code": "ESYNTAXUC", "details": "Unexpected character."}, **e.serialize())


def syntax_ueof(e):
    return dict({"code": "ESYNTAXUEOF", "details": "Unexpected end of input."}, **e.serialize())


def syntax_nested(e):
    return {"code": "ESYNTAXNESTED", "details": "Nested descriptions were encountered during parsing."}


def position_syntax(details, e):
    return dict({"code": "EPOSITIONSYNTAX", "details": details}, **e.serialize())


def position_invalid():
    # Used in the position converter.
    return {"code": "EPOSITIONINVALID", "details": "Position is invalid; it must be a string."}


def no_inputs():
    # Used in the position converter.
    return {"code": "ENOINPUTS"}


def no_inputs_other():
    # Used in the position converter.
    return {"code": "ENOINPUTSOTHER"}


def no_to_selector(reference_id, selector_id):
    # Used in the position converter.
    return {
        "code": "ENOTOSELECTOR",
        "details": f"No selector `{selector_id}` was found in `{reference_id}`.",
    }


def splice_site(path):
    return {
        "code": "ESPLICESITE",
        "details": "Variant affects one or more splice sites.",
        "paths": [path],
    }


def invalid_input(value, valid_options=None):
    # Used in the algebra compare input.
    output = {"code": "EINVALIDINPUT", "details": f"`{value}` is invalid."}
    if valid_options:
        output["options"] = valid_options
    return output


def sequence_length(seq, len_max):
    # Used in the mapper.
    return {
        "code": "ESEQUENCELENGTH",
        "details": f"Sequence mapping is not supported; sequence length {len(seq)} exceeds the maximum supported "
                   f"({len_max})."
    }


def slice_option(slice_to):
    # Used in the mapper.
    return {
        "code": "ESLICEOPTION",
        "details": f"Slicing to `{slice_to}` is not supported; only `gene` and `transcript` are supported."
    }


def location_slice(location):
    # Used in the mapper.
    return {
        "code": "ELOCATIONSLICE",
        "details": f"`{location_to_description(location)}` cannot be sliced for sequence mapping."
    }


def lengths_difference(length, accepted_difference):
    # Used in the mapper.
    return {
        "code": "ELENGTHSDIFFERENCE",
        "details": f"Sequence mapping is not allowed; length difference `{length}` bases exceeds the configured limit "
                   f"(`{accepted_difference}` bases)."
    }


def cds_slices(selector_id, exception_message):
    # NM_004152.3:c.205_685del
    return {
        "code": "ECDSSLICES",
        "details": f"`{selector_id}` has non-consecutive CDS slices; annotation exception: {exception_message}."
    }


def no_cds(reference_id, selector_id, path):
    # NG_009930.1(NM_001099625.2):c.1010
    return {
        "code": "ENOCDS",
        "details": f"`{reference_id}` has no CDS annotation for transcript selector `{selector_id}`.",
        "paths": [path],
    }
