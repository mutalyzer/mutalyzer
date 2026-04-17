from mutalyzer_crossmapper.crossmapper import Coding, Genomic, NonCoding
from mutalyzer_hgvs_parser import to_model
from mutalyzer_hgvs_parser.exceptions import UnexpectedCharacter, UnexpectedEnd
import mutalyzer.errors as errors
from mutalyzer.reference import retrieve_reference, get_internal_selector_model


REGIONS = {"", "u", "d", "*", "-"}

def validate_position_model(position_m_type: str, position_model: dict) -> None:
    """Validate the position model based on the coordinate system."""

    # TODO: check for positive value
    if position_m_type == "dna":
        # ignore offset, region and position_in_codon
        if not isinstance(position_model.get("position"), int):
            raise ValueError("Position must be an integer for genomic coordinate system.")

    elif position_m_type == "rna":
        # ignore position_in_codon
        if not isinstance(position_model.get("position"), int):
            raise ValueError("Position must be an integer for transcript coordinate system.")
        if not isinstance(position_model.get("offset"), int):
            raise ValueError("Offset must be an integer for transcript coordinate system.")
        if position_model.get("region") not in REGIONS:
            raise ValueError(f"Region must be one of {REGIONS} for transcript coordinate system.")

    elif position_m_type == "protein":
        if not isinstance(position_model.get("position"), int):
            raise ValueError("Position must be an integer for protein coordinate system.")
        if position_model.get("position_in_codon") not in (1, 2, 3):
            raise ValueError("Position in codon must be 1, 2, or 3.")
        if not isinstance(position_model.get("offset"), int):
            raise ValueError("Offset must be an integer for protein coordinate system.")
        if position_model.get("region") != "":
            raise ValueError(f"Region must be empty for protein coordinate system.")
    else:
        raise ValueError("Invalid HGVS position model type.")


def check_ref_length(reference_m: dict, coordinate: int) -> None:
    """Check if the coordinate exceeds the reference length."""
    reference_length = len(reference_m.get("sequence", {}).get("seq", ""))
    if coordinate > reference_length-1:
        raise ValueError(f"Coordinate {coordinate} exceeds reference length of {reference_length}.")
    if coordinate < 0:
        raise ValueError(f"Coordinate {coordinate} is less than 0.")


def validate_model(model: dict, r_id: str, s_id: str="") -> None:
    """Validate reference model or selector model, raise an error if not."""
    if not model:
        if s_id:
            raise ValueError(f"Selector {s_id} is not a valid selector in reference {r_id}.")
        raise ValueError(f"Reference {r_id} not retrieved.")


def retrieve_valid_reference_model(reference_id: str) -> dict:
    """Retrieve a reference model and validate it."""
    reference_m = retrieve_reference(reference_id)[0]
    validate_model(reference_m, reference_id)
    return reference_m


def resolve_selector_id(reference_id: str, selector_id: str) -> str:
    """Resolve an empty selector id to the reference id."""
    return reference_id if selector_id == "" else selector_id


def get_valid_selector_model(reference_m: dict, selector_id: str, reference_id: str, selector_type: str="") -> dict:
    """Retrieve and validate an internal selector model."""
    internal_selector_m = get_internal_selector_model(reference_m["annotations"], selector_id)
    validate_model(internal_selector_m, reference_id, selector_id)
    if selector_type:
        validate_selector_type(internal_selector_m, selector_type, selector_id)
    return internal_selector_m


def make_coding_crossmap(selector_m: dict) -> Coding:
    """Create a Coding crossmap from a selector model."""
    return Coding(selector_m["exon"], selector_m["cds"][0], selector_m["inverted"])


def validate_selector_type(model: dict, selector_type: str, selector_id: str) -> None:
    """Validate if the selector type is correct, raise an error if not."""
    if selector_type == "rna" and "rna" not in model.get("type", "").lower():
        raise ValueError(f"Selector {selector_id} is not a transcript.")
    if selector_type == "protein" and "cds" not in model.get("type", "").lower():
        raise ValueError(f"Selector {selector_id} is not a protein.")
    # support for gene as a selector?


def validate_coding_transcript(model: dict, trancript_id: str) -> None:
    """Validate if the selector is a coding transcript, raise an error if not."""
    if not model.get("cds"):
        raise ValueError(f"Selector {trancript_id} is not a coding transcript.")


def check_intron_exon_boundary(crossmap: Coding, exon_list: list, position_m: dict) -> None:
    """Check if the position is at the intron-exon boundary, raise an error if not."""
    exon_start = [exon[0] for exon in exon_list]
    exon_end = [exon[1] for exon in exon_list]
    position = position_m['position']
    offset = position_m['offset']

    exon_pos_model = {'position': position, 'offset': 0, 'region': position_m['region']}

    coordinate = crossmap.coding_to_coordinate(exon_pos_model)
    # TODO: check for reverse strand
    if coordinate in exon_start:
        if offset < 0:
            return
        else:
            raise ValueError(f"{position} is at the exon start.")
    if coordinate + 1 in exon_end:
        if offset > 0:
            return
        else:
            raise ValueError(f"{position} is at the exon end.")
    raise ValueError(f"Position {position} is not on the intron-exon boundary.")


def coordinate_to_genomic(coordinate: int) -> dict:
    """Convert a coordinate to HGVS genomic position.

    Args:
        coordinate (int): coordinate position to convert.

    Returns:
        dict: The converted HGVS genomic position model.
    """
    crossmap = Genomic()
    return crossmap.coordinate_to_genomic(coordinate)


def coordinate_to_coding(reference_id: str, coordinate: int, transcript_id: str="") -> dict:
    """Convert a coordinate to HGVS coding position model on a genomic reference sequence.

    Args:
        reference_id (str): ID of the reference sequence.
        coordinate (int): coordinate to convert.
        transcript_id (str): ID of the transcript sequence.

    Returns:
        dict: The converted HGVS coding position model.
    """
    reference_m = retrieve_valid_reference_model(reference_id)
    check_ref_length(reference_m, coordinate)
    if transcript_id:
        selector_id = transcript_id
    else:
        selector_id = reference_id
    internal_selector_m = get_valid_selector_model(reference_m, selector_id, reference_id, selector_type="rna")
    validate_coding_transcript(internal_selector_m, selector_id)
    return make_coding_crossmap(internal_selector_m).coordinate_to_coding(coordinate)


def coordinate_to_noncoding(reference_id: str, coordinate: int, transcript_id: str="") -> dict:
    """Convert a coordinate to HGVS non-coding position model on a genomic reference sequence.

    Args:
        reference_id (str): ID of the reference sequence.
        coordinate (int): coordinate to convert.
        transcript_id (str): ID of the non-coding transcript sequence.

    Returns:
        dict: The converted position model in the HGVS non-coding coordinate system.
    """
    reference_m = retrieve_valid_reference_model(reference_id)
    check_ref_length(reference_m, coordinate)
    if transcript_id:
        selector_id = transcript_id
    else:
        selector_id = reference_id

    internal_selector_model = get_valid_selector_model(reference_m, selector_id, reference_id, selector_type="rna")

    # Allow for coding transcript and return noncoding position model
    crossmap = NonCoding(internal_selector_model["exon"], internal_selector_model["inverted"])
    return crossmap.coordinate_to_noncoding(coordinate)


def coordinate_to_protein(coordinate: int) -> dict:
    """Convert a coordinate to HGVS protein position model on a protein sequence.

    Args:
        coordinate (int): coordinate to convert.

    Returns
        dict: The converted position model in the HGVS protein coordinate system.
    """
    position_m = coordinate_to_genomic(coordinate)
    # QUESTION: position_in_codon value
    position_m['position_in_codon'] = 1
    position_m['offset'] = 0
    position_m['region'] = ""

    return position_m


def coordinate_to_reference_protein(reference_id: str, coordinate: int, protein_id: str) -> dict:
    """Convert a coordinate to protein coordinate system on a genomic reference sequence.

    Args:
        reference_id (str):  ID of the reference sequence.
        coordinate (int): coordinate to convert.
        protein_id (str):  ID of the protein sequence.

    Returns:
        dict: The converted position model in the HGVS protein coordinate system.
    """
    reference_m = retrieve_valid_reference_model(reference_id)
    check_ref_length(reference_m, coordinate)
    internal_protein_m = get_valid_selector_model(reference_m, protein_id, reference_id, selector_type="protein")
    return make_coding_crossmap(internal_protein_m).coordinate_to_protein(coordinate)


def genomic_to_coordinate(position_model: dict) -> int:
    """Convert a HGVS genomic position model to coordinate.

    Args:
        position_model (dict): The HGVS genomic position model to convert.

    Returns:
        int: The converted coordinate.
    """
    validate_position_model("dna", position_model)
    crossmap = Genomic()
    return crossmap.genomic_to_coordinate(position_model)


def genomic_coding_to_coordinate(reference_id: str, transcript_id: str, position_m: dict) -> int:
    """Convert a HGVS coding position model on a genomic reference sequence to coordinate.

    Args:
        reference_id (str): ID of the reference sequence.
        transcript_id (str): ID of the transcript sequence.
        position_model (dict): The HGVS coding position model to convert.

    Returns:
        int: The converted coordinate.
    """
    validate_position_model("rna", position_m)
    reference_m = retrieve_valid_reference_model(reference_id)
    internal_selector_m = get_valid_selector_model(reference_m, transcript_id, reference_id, selector_type="rna")
    validate_coding_transcript(internal_selector_m, transcript_id)
    crossmap = Coding(internal_selector_m["exon"], internal_selector_m["cds"][0], internal_selector_m["inverted"])
    coordinate = crossmap.coding_to_coordinate(position_m)

    if position_m.get("offset") != 0:
        check_intron_exon_boundary(crossmap, internal_selector_m["exon"], position_m)
    check_ref_length(reference_m, coordinate)
    return coordinate


def transcript_coding_to_coordinate(transcript_id: str, position_m: dict) -> int:
    """Convert a HGVS coding position model on a transcript to coordinate.

    Args:
        transcript_id (str): ID of the transcript sequence.
        position_model (dict): The HGVS coding position model to convert.

    Returns:
        int: The converted coordinate.
    """
    validate_position_model("rna", position_m)
    reference_m = retrieve_valid_reference_model(transcript_id)
    internal_selector_m = get_valid_selector_model(reference_m, transcript_id, transcript_id, selector_type="rna")
    validate_coding_transcript(internal_selector_m, transcript_id)
    crossmap = Coding(internal_selector_m["exon"], internal_selector_m["cds"][0], internal_selector_m["inverted"])
    coordinate = crossmap.coding_to_coordinate(position_m)
    check_ref_length(reference_m, coordinate)
    return coordinate


def genomic_noncoding_to_coordinate(reference_id: str, transcript_id: str, position_m: dict) -> int:
    """Convert a HGVS non-coding position model on a genomic reference sequence to coordinate.

    Args:
        reference_id (str): ID of the reference sequence.
        transcript_id (str): ID of the non-coding transcript sequence.
        position_model (dict): The HGVS non-coding position model to convert.

    Returns:
        int: The converted coordinate.
    """
    validate_position_model("rna", position_m)
    reference_m = retrieve_valid_reference_model(reference_id)
    internal_selector_m = get_valid_selector_model(reference_m, transcript_id, reference_id, selector_type="rna")
    crossmap = NonCoding(internal_selector_m["exon"], internal_selector_m["inverted"])

    if position_m.get("offset") != 0:
        check_intron_exon_boundary(crossmap, internal_selector_m["exon"], position_m)
    coordinate = crossmap.noncoding_to_coordinate(position_m)
    check_ref_length(reference_m, coordinate)
    return coordinate


def transcript_noncoding_to_coordinate(transcript_id: str, position_m: dict) -> int:
    """Convert a HGVS non-coding position model on a transcript to coordinate.

    Args:
        transcript_id (str): ID of the non-coding transcript sequence.
        position_model (dict): The HGVS non-coding position model to convert.

    Returns:
        int: The converted coordinate.
    """
    validate_position_model("rna", position_m)
    reference_m = retrieve_valid_reference_model(transcript_id)
    internal_selector_m = get_valid_selector_model(reference_m, transcript_id, transcript_id, selector_type="rna")
    crossmap = NonCoding(internal_selector_m["exon"], internal_selector_m["inverted"])
    coordinate = crossmap.noncoding_to_coordinate(position_m)
    check_ref_length(reference_m, coordinate)
    return coordinate


def protein_to_coordinate(position_m: dict) -> int:
    """Convert a HGVS protein position model on a protein sequence to coordinate.

    Args:
        protein_id (str): ID of the protein sequence.
        position_model (dict): The HGVS protein position model to convert.

    Returns:
        int: The converted coordinate.
    """
    validate_position_model("protein", position_m)
    return genomic_to_coordinate(position_m)


def transcript_protein_to_coordinate(reference_id: str, protein_id: str, position_m: dict) -> int:
    """Convert a HGVS protein position model on a genomic reference sequence to coordinate.

    Args:
        reference_id (str): ID of the reference sequence.
        protein_id (str): ID of the protein sequence.
        position_model (dict): The HGVS protein position model to convert.

    Returns:
        int: The converted coordinate.
    """
    validate_position_model("protein", position_m)
    reference_m = retrieve_valid_reference_model(reference_id)
    internal_protein_m = get_valid_selector_model(reference_m, protein_id, reference_id, selector_type="protein")
    crossmap = Coding(internal_protein_m["exon"], internal_protein_m["cds"][0], internal_protein_m["inverted"])
    coordinate = crossmap.noncoding_to_coordinate(position_m)
    check_ref_length(reference_m, coordinate)
    return make_coding_crossmap(internal_protein_m).protein_to_coordinate(position_m)


def genomic_protein_to_coordinate(reference_id: str, protein_id: str, position_m: dict) -> int:
    """Convert a HGVS protein position model on a genomic reference sequence to coordinate.

    Args:
        reference_id (str): ID of the reference sequence.
        protein_id (str): ID of the protein sequence.
        position_model (dict): The HGVS protein position model to convert.

    Returns:
        int: The converted coordinate.
    """
    validate_position_model("protein", position_m)
    reference_m = retrieve_valid_reference_model(reference_id)
    internal_protein_m = get_valid_selector_model(reference_m, protein_id, reference_id, selector_type="protein")
    crossmap = Coding(internal_protein_m["exon"], internal_protein_m["cds"][0], internal_protein_m["inverted"])
    coordinate = crossmap.noncoding_to_coordinate(position_m)
    check_ref_length(reference_m, coordinate)
    return make_coding_crossmap(internal_protein_m).protein_to_coordinate(position_m)
