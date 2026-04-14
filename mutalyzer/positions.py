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



def check_ref_length(reference_model: dict, position: int) -> None:
    """Check if the position exceeds the reference length."""
    reference_length = len(reference_model.get("sequence", {}).get("seq", ""))
    if position > reference_length:
        raise ValueError(f"Position {position} exceeds reference length of {reference_length}.")

def check_ref(reference_model: dict) -> None:
    """Check if the reference model is valid."""
    if not reference_model:
        raise ValueError("Reference not found.")



def coordinate_to_genomic_coding(reference_id: str, position: int, selector_id: str) -> dict:
    """Convert a standard coordinate to a coding position system on genomic reference sequence.

    Args:
        reference_id (str): The NCBI ID of the reference sequence.
        position (int): The coordinate position to convert.
        selector_id (str): The NCBI ID of the selector sequence.

    Returns
        dict: The converted position model in the HGVS coding coordinate system.
    """
    reference_model = retrieve_reference(reference_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, position)

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

    position_m = crossmap.coordinate_to_coding(position)
    return position_m


def coordinate_to_transcript_coding(transcript_id: str, position: int) -> dict:
    """Convert a standard coordinate to a coding position system on a transcript sequence.

    Args:
        transcript_id (str): The NCBI ID of the transcript sequence.
        position (int): The coordinate position to convert.

    Returns:
        dict: The converted position model in the HGVS coding coordinate system.
    """
    reference_model = retrieve_reference(transcript_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, position)

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

    position_m = crossmap.coordinate_to_coding(position)
    return position_m


def coordinate_to_genomic(reference_id: str, position: int) -> dict:
    """Convert a standard coordinate to genomic coordinate system.

    Args:
        reference_id (str): The NCBI ID of the reference sequence.
        position (int): The coordinate position to convert.

    Returns
        dict: The converted position model in the HGVS genomic coordinate system.
    """
    reference_model = retrieve_reference(reference_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, position)

    crossmap = Genomic()
    position_m = crossmap.coordinate_to_genomic(position)

    return position_m


def coordinate_to_genomic_noncoding(reference_id: str, position: int, selector_id: str) -> dict:
    """Convert a standard coordinate to non-coding position system on a genomic reference sequence.

    Args:
        reference_id (str): The NCBI ID of the reference sequence.
        position (int): The coordinate position to convert.
        selector_id (str): The NCBI ID of the non-codingselector sequence.

    Returns
        dict: The converted position model in the HGVS non-coding coordinate system.
    """
    reference_model = retrieve_reference(reference_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, position)

    selector_model = get_internal_selector_model(reference_model["annotations"], selector_id=selector_id)
    if not is_selector_in_reference(selector_id, reference_model):
        raise ValueError(f"Selector {selector_id} not found in reference {reference_id}.")


    # Allow for coding transcript and return noncoding position model
    crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
    position_m = crossmap.coordinate_to_noncoding(position)

    return position_m


def coordinate_to_transcript_noncoding(transcript_id: str, position: int) -> dict:
    """Convert a standard coordinate to non-coding position system on a transcript.

    Args:
        transcript_id (str): The NCBI ID of the non-coding transcript sequence.
        position (int): The coordinate position to convert.

    Returns
        dict: The converted position model in the HGVS non-coding coordinate system.
    """
    reference_model = retrieve_reference(transcript_id)[0]
    check_ref(reference_model)
    check_ref_length(reference_model, position)

    selector_model = get_internal_selector_model(reference_model["annotations"], selector_id=transcript_id)
    if not is_selector_in_reference(transcript_id, reference_model):
        raise ValueError(f"Selector {transcript_id} not found in reference {transcript_id}.")

    # Allow for coding transcript and return noncoding position model
    crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
    position_m = crossmap.coordinate_to_noncoding(position)

    return position_m
