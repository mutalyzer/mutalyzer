from mutalyzer_crossmapper.crossmapper import Coding, Genomic, NonCoding
from mutalyzer.position_converter import position_convert
# from .converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer_hgvs_parser import to_model
from mutalyzer_hgvs_parser.exceptions import UnexpectedCharacter, UnexpectedEnd
import mutalyzer.errors as errors
from mutalyzer.reference import retrieve_reference, get_internal_selector_model


ALLOWED_REGIONS = {"", "u", "d", "*", "-"}

def validate_position_model(coordinate_system: str, position_model: dict) -> None:
    """Validate the position model based on the coordinate system."""
    if coordinate_system == "genomic":
        if not isinstance(position_model.get("position"), int):
            raise ValueError("Position must be an integer for genomic coordinate system.")

    elif coordinate_system == "transcript":
        if not isinstance(position_model.get("position"), int):
            raise ValueError("Position must be an integer for transcript coordinate system.")
        if not isinstance(position_model.get("offset"), int):
            raise ValueError("Offset must be an integer for transcript coordinate system.")
        if position_model.get("region") not in ALLOWED_REGIONS:
            raise ValueError(f"Region must be one of {ALLOWED_REGIONS} for transcript coordinate system.")

    elif coordinate_system == "protein":
        if not isinstance(position_model.get("position"), int):
            raise ValueError("Position must be an integer for protein coordinate system.")
        if position_model.get("position_in_codon") not in (1, 2, 3):
            raise ValueError("Position in codon must be 1, 2, or 3.")
        if not isinstance(position_model.get("offset"), int):
            raise ValueError("Offset must be an integer for protein coordinate system.")
        if position_model.get("region") not in ALLOWED_REGIONS:
            raise ValueError(f"Region must be one of {ALLOWED_REGIONS} for protein coordinate system.")
    else:
        raise ValueError("Invalid coordinate system.")


def check_ref_length(reference_m: dict, coordinate: int) -> None:
    """Check if the coordinate exceeds the reference length."""
    reference_length = len(reference_m.get("sequence", {}).get("seq", ""))
    if coordinate > reference_length-1:
        raise ValueError(f"Coordinate {coordinate} exceeds reference length of {reference_length}.")


def validate_model(model: dict, r_id: str, s_id: str="") -> None:
    """Validate reference model or selector model, raise an error if not."""
    if not model:
        if s_id:
            raise ValueError(f"{s_id} not found in reference {r_id}.")
        raise ValueError(f"{r_id} not retrieved.")


def validate_selector_type(model: dict, bio_type: str, selector_id: str) -> None:
    """Validate if the selector type is correct, raise an error if not."""
    if bio_type == "rna" and "rna" not in model.get("type", "").lower():
        raise ValueError(f"{selector_id} is not a transcript.")
    if bio_type == "protein" and "cds" not in model.get("type", "").lower():
        raise ValueError(f"{selector_id} is not a protein.")
    # support for gene as a selector?


def validate_coding_transcript(model: dict, trancript_id: str) -> None:
    """Validate if the selector is a coding transcript, raise an error if not."""
    if not model.get("cds"):
        raise ValueError(f"{trancript_id} is not a coding transcript.")


def coordinate_to_genomic(coodinate: int) -> dict:
    """Convert a coordinate to HGVS genomic position.

    Args:
        coordinate (int): coordinate position to convert.

    Returnscoding position system on a transcript sequence
        dict: The converted HGVS genomic position model.
    """
    crossmap = Genomic()
    return crossmap.coordinate_to_genomic(coodinate)


def coordinate_to_genomic_coding(reference_id: str, coordinate: int, transcript_id: str) -> dict:
    """Convert a coordinate to HGVS coding position model on a genomic reference sequence.

    Args:
        reference_id (str): ID of the reference sequence.
        coordinate (int): coordinate to convert.
        transcript_id (str): ID of the transcript sequence.

    Returns
        dict: The converted HGVS coding position model.
    """
    reference_m = retrieve_reference(reference_id)[0]
    validate_model(reference_m, reference_id)
    check_ref_length(reference_m, coordinate)

    internal_selector_m = get_internal_selector_model(reference_m["annotations"], selector_id=transcript_id)
    validate_model(internal_selector_m, reference_id, transcript_id)
    validate_selector_type(internal_selector_m, "rna", transcript_id)
    validate_coding_transcript(internal_selector_m, transcript_id)

    crossmap = Coding(
        internal_selector_m["exon"],
        internal_selector_m["cds"][0],
        internal_selector_m["inverted"]
    )

    return crossmap.coordinate_to_coding(coordinate)


def coordinate_to_transcript_coding(transcript_id: str, coordinate: int) -> dict:
    """Convert a coordinate to a HGVS coding position model on a transcript.

    Args:
        transcript_id (str): ID of the transcript sequence.
        coordinate (int): coordinate to convert.

    Returns:
        dict: The converted HGVS coding position model.
    """
    reference_m = retrieve_reference(transcript_id)[0]
    validate_model(reference_m, transcript_id)
    check_ref_length(reference_m, coordinate)

    internal_selector_m = get_internal_selector_model(reference_m["annotations"], transcript_id)
    validate_selector_type(internal_selector_m, "rna", transcript_id)
    validate_coding_transcript(internal_selector_m, transcript_id)

    crossmap = Coding(
        internal_selector_m["exon"],
        internal_selector_m["cds"][0],
        internal_selector_m["inverted"]
    )

    return crossmap.coordinate_to_coding(coordinate)


def coordinate_to_genomic_noncoding(reference_id: str, coordinate: int, transcript_id: str) -> dict:
    """Convert a coordinate to HGVS non-coding position model on a genomic reference sequence.

    Args:
        reference_id (str): ID of the reference sequence.
        coordinate (int): coordinate to convert.
        transcript_id (str): ID of the non-coding transcript sequence.

    Returns
        dict: The converted position model in the HGVS non-coding coordinate system.
    """
    reference_m = retrieve_reference(reference_id)[0]
    validate_model(reference_m, reference_id)
    check_ref_length(reference_m, coordinate)

    internal_selector_model = get_internal_selector_model(reference_m["annotations"], transcript_id)
    validate_model(internal_selector_model, reference_id, transcript_id)
    validate_selector_type(internal_selector_model, "rna", transcript_id)

    # Allow for coding transcript and return noncoding position model
    crossmap = NonCoding(internal_selector_model["exon"], internal_selector_model["inverted"])

    return crossmap.coordinate_to_noncoding(coordinate)


def coordinate_to_transcript_noncoding(transcript_id: str, coordinate: int) -> dict:
    """Convert a coordinate to HGVS non-coding position model on a transcript.

    Args:
        transcript_id (str): ID of the non-coding transcript sequence.
        coordinate (int): coordinate to convert.

    Returns
        dict: The converted position model in the HGVS non-coding coordinate system.
    """
    transcript_m = retrieve_reference(transcript_id)[0]
    validate_model(transcript_m, transcript_id)
    check_ref_length(transcript_m, coordinate)

    internal_selector_m = get_internal_selector_model(transcript_m["annotations"], transcript_id)
    validate_selector_type(internal_selector_m, "rna", transcript_id)

    # Allow for coding transcript and return noncoding position model
    crossmap = NonCoding(internal_selector_m["exon"], internal_selector_m["inverted"])

    return crossmap.coordinate_to_noncoding(coordinate)


def coordinate_to_protein(coordinate: int) -> dict:
    """Convert a coordinate to HGVS protein position model on a protein sequence.

    Args:
        coordinate (int): coordinate to convert.

    Returns
        dict: The converted position model in the HGVS protein coordinate system.
    """
    crossmap = Genomic()
    position_m = crossmap.coordinate_to_genomic(coordinate)
    position_m["position_in_codon"] = 1
    position_m["offset"] = 0
    position_m["region"] = ""

    return position_m


def coordinate_to_reference_protein(reference_id: str, coordinate: int, protein_id: str) -> dict:
    """Convert a coordinate to protein coordinate system on a genomic reference sequence.

    Args:
        reference_id (str):  ID of the reference sequence.
        coordinate (int): coordinate to convert.
        protein_id (str):  ID of the protein sequence.

    Returns
        dict: The converted position model in the HGVS protein coordinate system.
    """
    reference_m = retrieve_reference(reference_id)[0]
    validate_model(reference_m, reference_id)
    check_ref_length(reference_m, coordinate)

    internal_protein_m = get_internal_selector_model(reference_m["annotations"], protein_id)
    validate_model(internal_protein_m, reference_id, protein_id)
    validate_selector_type(internal_protein_m, "protein", protein_id)

    crossmap = Coding(
        internal_protein_m["exon"],
        internal_protein_m["cds"][0],
        internal_protein_m["inverted"]
    )

    return crossmap.coordinate_to_protein(coordinate)
