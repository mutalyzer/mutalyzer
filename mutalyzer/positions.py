import argparse
import pprint

from mutalyzer_crossmapper.crossmapper import Coding, Genomic, NonCoding
from mutalyzer.position_converter import position_convert
# from .converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer_hgvs_parser import to_model
from mutalyzer_hgvs_parser.exceptions import UnexpectedCharacter, UnexpectedEnd
import mutalyzer.errors as errors
from mutalyzer.reference import retrieve_reference, get_internal_selector_model, is_selector_in_reference
from mutalyzer.hgvs_position_model import HGVSPositionModel


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

def check_ref_length(reference_model: dict, coordinate: int) -> None:
    """Check if the coordinate exceeds the reference length."""
    reference_length = len(reference_model.get("sequence", {}).get("seq", ""))
    if coordinate > reference_length-1:
        raise ValueError(f"Coordinate {coordinate} exceeds reference length of {reference_length}.")

def check_ref(reference_model: dict) -> None:
    """Check if the reference model is valid."""
    if not reference_model:
        raise ValueError("Reference not found.")



def coordinate_to_genomic_coding(reference_id: str, coordinate: int, selector_id: str) -> dict:
    """Convert a standard coordinate to a coding position system on genomic reference sequence.

    Args:
        reference_id (str): The NCBI ID of the reference sequence.
        coordinate (int): The coordinate to convert.
        selector_id (str): The NCBI ID of the selector sequence.

    Returns
        dict: The converted position model in the HGVS coding coordinate system.
    """
    reference_model = retrieve_reference(reference_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, coordinate)

    selector_model = get_internal_selector_model(reference_model["annotations"], selector_id=selector_id)
    if not is_selector_in_reference(selector_id, reference_model):
        raise ValueError(f"Selector {selector_id} not found in reference {reference_id}.")
    if not selector_model.get("cds"):
        raise ValueError(f"Transcript {selector_id} is not a coding transcript, cannot convert to coding coordinate.")

    crossmap = Coding(
        selector_model["exon"],
        selector_model["cds"][0],
        selector_model["inverted"]
    )

    position_m = crossmap.coordinate_to_coding(coordinate)
    return position_m


def genomic_coding_to_coordinate(reference_id: str, position_model: dict, selector_id: str) -> int:
    """Convert a coding position system on genomic reference sequence to standard coordinate.

    Args:
        reference_id (str): The NCBI ID of the reference sequence.
        position_model (dict): The position model in the HGVS coding coordinate system to convert.
        selector_id (str): The NCBI ID of the selector sequence.

    Returns
        int: Standard 0-based coordinate.
    """

    position = position_model.get("position")
    if position is None:
        raise ValueError("Position is required in position model.")

    reference_model = retrieve_reference(reference_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, abs(position))

    selector_model = get_internal_selector_model(reference_model["annotations"], selector_id=selector_id)
    if not is_selector_in_reference(selector_id, reference_model):
        raise ValueError(f"Selector {selector_id} not found in reference {reference_id}.")
    if not selector_model.get("cds"):
        raise ValueError(f"Transcript {selector_id} is not a coding transcript, cannot convert from coding coordinate.")

    crossmap = Coding(
        selector_model["exon"],
        selector_model["cds"][0],
        selector_model["inverted"]
    )

    coordinate = crossmap.coding_to_coordinate(position_model)
    return coordinate


def coordinate_to_transcript_coding(transcript_id: str, coordinate: int) -> dict:
    """Convert a standard coordinate to a coding position system on a transcript sequence.

    Args:
        transcript_id (str): The NCBI ID of the transcript sequence.
        coordinate (int): The coordinate to convert.

    Returns:
        dict: The converted position model in the HGVS coding coordinate system.
    """
    reference_model = retrieve_reference(transcript_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, coordinate)

    selector_model = get_internal_selector_model(reference_model["annotations"], selector_id=transcript_id)
    if not is_selector_in_reference(transcript_id, reference_model):
        raise ValueError(f"Transcript {transcript_id} not found in reference {transcript_id}.")
    if not selector_model.get("cds"):
        raise ValueError(f"Transcript {transcript_id} is not a coding transcript, cannot convert to coding coordinate.")

    crossmap = Coding(
        selector_model["exon"],
        selector_model["cds"][0],
        selector_model["inverted"]
    )

    position_m = crossmap.coordinate_to_coding(coordinate)
    return position_m


def transcript_coding_to_coordinate(transcript_id: str, position_model: dict) -> int:
    """Convert a coding position system on a transcript sequence to standard coordinate.

    Args:
        transcript_id (str): The NCBI ID of the transcript sequence.
        position_model (dict): The position model in the HGVS coding coordinate system to convert.

    Returns:
        int: Standard 0-based coordinate.
    """
    position = position_model.get("position")
    if position is None:
        raise ValueError("Position is required in position model.")

    reference_model = retrieve_reference(transcript_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, abs(position))
    # check if position model is valid

def coordinate_to_genomic(coodinate: int) -> dict:
    """Convert a standard coordinate to genomic coordinate system.

    Args:
        coordinate (int): The coordinate position to convert.

    Returns
        dict: The converted position model in the HGVS genomic coordinate system.
    """
    crossmap = Genomic()
    position_m = crossmap.coordinate_to_genomic(coodinate)

    return position_m


def coordinate_to_genomic_noncoding(reference_id: str, coordinate: int, selector_id: str) -> dict:
    """Convert a standard coordinate to non-coding position system on a genomic reference sequence.

    Args:
        reference_id (str): The NCBI ID of the reference sequence.
        coordinate (int): The coordinate to convert.
        selector_id (str): The NCBI ID of the non-codingselector sequence.

    Returns
        dict: The converted position model in the HGVS non-coding coordinate system.
    """
    reference_model = retrieve_reference(reference_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, coordinate)

    selector_model = get_internal_selector_model(reference_model["annotations"], selector_id=selector_id)
    if not is_selector_in_reference(selector_id, reference_model):
        raise ValueError(f"Selector {selector_id} not found in reference {reference_id}.")

    # Allow for coding transcript and return noncoding position model
    crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
    position_m = crossmap.coordinate_to_noncoding(coordinate)

    return position_m


def coordinate_to_transcript_noncoding(transcript_id: str, coordinate: int) -> dict:
    """Convert a standard coordinate to non-coding position system on a transcript.

    Args:
        transcript_id (str): The NCBI ID of the non-coding transcript sequence.
        coordinate (int): The coordinate to convert.

    Returns
        dict: The converted position model in the HGVS non-coding coordinate system.
    """
    reference_model = retrieve_reference(transcript_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, coordinate)

    selector_model = get_internal_selector_model(reference_model["annotations"], selector_id=transcript_id)
    if not is_selector_in_reference(transcript_id, reference_model):
        raise ValueError(f"Selector {transcript_id} not found in reference {transcript_id}.")

    # Allow for coding transcript and return noncoding position model
    crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
    position_m = crossmap.coordinate_to_noncoding(coordinate)

    return position_m


def coordinate_to_protein(coordinate: int) -> dict:
    """Convert a standard coordinate to protein coordinate system.

    Args:
        coordinate (int): The coordinate to convert.

    Returns
        dict: The converted position model in the HGVS protein coordinate system.
    """
    crossmap = Genomic()
    position_m = crossmap.coordinate_to_genomic(coordinate)
    position_m["position_in_codon"] = 0
    position_m["offset"] = 0
    position_m["region"] = ""

    return position_m


def coordinate_to_transcript_protein(transcript_id: str, coordinate: int, protein_id: str) -> dict:
    """Convert a standard coordinate to protein coordinate system on a transcript.

    Args:
        transcript_id (str): The NCBI ID of the transcript sequence.
        coordinate (int): The coordinate to convert.
        protein_id (str): The NCBI ID of the protein sequence.

    Returns
        dict: The converted position model in the HGVS protein coordinate system.
    """
    reference_model = retrieve_reference(transcript_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, coordinate)
    transcript_model = get_internal_selector_model(reference_model["annotations"], selector_id=transcript_id)
    if not transcript_model.get("cds"):
        raise ValueError(f"Transcript {transcript_id} is not a coding transcript, cannot convert to protein coordinate.")

    if not is_selector_in_reference(protein_id, reference_model):
        raise ValueError(f"Transcript {protein_id} not found in reference {transcript_id}.")

    crossmap = Coding(
        transcript_model["exon"],
        transcript_model["cds"][0],
        transcript_model["inverted"]
    )
    position_m = crossmap.coordinate_to_protein(coordinate)
    return position_m



def coordinate_to_genomic_protein(reference_id: str, coordinate: int, protein_id: str) -> dict:
    """Convert a standard coordinate to protein coordinate system on a genomic reference sequence.

    Args:
        reference_id (str): The NCBI ID of the reference sequence.
        coordinate (int): The coordinate to convert.
        protein_id (str): The NCBI ID of the protein sequence.

    Returns
        dict: The converted position model in the HGVS protein coordinate system.
    """
    reference_model = retrieve_reference(reference_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, coordinate)

    if not is_selector_in_reference(protein_id, reference_model):
        raise ValueError(f"Transcript {protein_id} not found in reference {reference_id}.")
    return reference_model