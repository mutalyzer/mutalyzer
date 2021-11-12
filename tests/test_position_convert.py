from mutalyzer_hgvs_parser import to_model

from mutalyzer.description_model import model_to_string
from mutalyzer.position_converter import position_convert

from .commons import code_in, get_codes, patch_retriever


def test_error_no_required_inputs_reference_id():
    p_c = position_convert(
        reference_id=None,
        position="100",
        from_selector_id=None,
        from_coordinate_system=None,
        to_coordinate_system=None,
        to_selector_id=None,
        include_overlapping=False,
    )
    assert get_codes(p_c["errors"]) == {"ENOINPUTS", "ENOINPUTSOTHER"}


def test_error_no_required_inputs_position():
    p_c = position_convert(
        reference_id="NG_012337.1",
        position=None,
        from_selector_id=None,
        from_coordinate_system=None,
        to_coordinate_system=None,
        to_selector_id=None,
        include_overlapping=False,
    )
    assert get_codes(p_c["errors"]) == {"ENOINPUTS", "ENOINPUTSOTHER"}


def test_error_no_required_inputs_other():
    p_c = position_convert(
        reference_id="NG_012337.1",
        position="100",
        from_selector_id=None,
        from_coordinate_system=None,
        to_coordinate_system=None,
        to_selector_id=None,
        include_overlapping=False,
    )
    assert p_c["errors"][0]["code"] == "ENOINPUTSOTHER"


def test_error_no_from_selector():
    p_c = position_convert(
        reference_id="NM_003002.4",
        from_selector_id="NM_",
        position="100",
    )
    assert p_c["errors"][0]["code"] == "ENOSELECTORFOUND"


def test_error_position_invalid():
    p_c = position_convert(
        reference_id="NM_003002.4",
        from_selector_id="NM_",
        position=100,
    )
    assert p_c["errors"][0]["code"] == "EPOSITIONINVALID"


def test_error_position_syntax():
    p_c = position_convert(
        reference_id="NG_012337.1",
        from_selector_id="NM_003002.2",
        from_coordinate_system="c",
        to_coordinate_system="g",
        position="100_",
    )
    assert p_c["errors"][0]["code"] == "EPOSITIONSYNTAX"


def test_info_from_selector_model_identified_from_selector_id():
    p_c = position_convert(
        reference_id="NG_012337.1",
        from_selector_id="NM_003002.2",
        position="100",
    )
    assert p_c["infos"][0]["code"] == "ICORRECTEDCOORDINATESYSTEM"


def test_info_from_selector_model_identifed_from_reference():
    p_c = position_convert(
        reference_id="NG_012337.1",
        position="100",
        to_selector_id="NM_003002.2",
    )
    assert p_c["infos"][0]["code"] == "ICORRECTEDCOORDINATESYSTEM"


def test_info_to_selector_model_identified_from_reference():
    p_c = position_convert(
        reference_id="NG_012337.1",
        from_selector_id="NM_003002.2",
        from_coordinate_system="c",
        position="100",
    )
    assert p_c["infos"][0]["code"] == "ITOSELECTORFROMREFERENCE"


def test_error_from_to_selectors_equal_g():
    p_c = position_convert(
        reference_id="NG_012337.1",
        from_coordinate_system="g",
        position="100",
        to_coordinate_system="g",
    )
    assert p_c["infos"][0]["code"] == "IFROMTOSELECTORSEQUAL"


def test_error_invalid_to_selector():
    p_c = position_convert(
        reference_id="NG_012337.1",
        position="100",
        to_selector_id="NM_",
    )
    assert p_c["errors"][0]["code"] == "ENOTOSELECTOR"


def test_outside_cds_segmented():
    p_c = position_convert(
        reference_id="NG_012337.1",
        from_selector_id=None,
        from_coordinate_system=None,
        position="27",
        to_coordinate_system=None,
        to_selector_id="NM_003002.2",
    )
    assert model_to_string(p_c["converted_model"]) == "NG_012337.1(NM_003002.2):c.-5035"


def test_outside_cds_description():
    p_c = position_convert(
        description="NG_012337.1:27",
        to_selector_id="NM_003002.2",
    )
    assert model_to_string(p_c["converted_model"]) == "NG_012337.1(NM_003002.2):c.-5035"


def test_outside_cds_description_model():
    p_c = position_convert(
        description_model=to_model("NG_012337.1:27"),
        to_selector_id="NM_003002.2",
    )
    assert model_to_string(p_c["converted_model"]) == "NG_012337.1(NM_003002.2):c.-5035"


def test_point_shift_from_positive_strand_to_negative_strand():
    """
    NM_012459.2 is on the negative strand on NG_012337.1
    """
    model = to_model("NG_012337.1:1005")
    model["variants"][0]["location"]["shift"] = 4
    p_c = position_convert(
        description_model=model,
        to_selector_id="NM_012459.2",
    )
    assert model_to_string(p_c["converted_model"]) == "NG_012337.1(NM_012459.2):c.*2448"


def test_point_shift_from_positive_strand_ro_positive_strand():
    model = to_model("NG_012337.1:1005")
    model["variants"][0]["location"]["shift"] = 4
    p_c = position_convert(
        description_model=model,
        to_selector_id="NM_003002.2",
    )
    assert model_to_string(p_c["converted_model"]) == "NG_012337.1(NM_003002.2):c.-4057"


def test_range_shift_from_positive_strand_to_negative_strand():
    """
    NM_012459.2 is on the negative strand on NG_012337.1
    """
    model = to_model("NG_012337.1:1004_1005")
    model["variants"][0]["location"]["start"]["shift"] = 3
    model["variants"][0]["location"]["end"]["shift"] = 3
    p_c = position_convert(
        description_model=model,
        to_selector_id="NM_012459.2",
    )
    assert (
        model_to_string(p_c["converted_model"])
        == "NG_012337.1(NM_012459.2):c.*2447_*2448"
    )
