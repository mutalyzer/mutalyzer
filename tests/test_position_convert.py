from mutalyzer_hgvs_parser import parse_description_to_model

from normalizer.description_model import model_to_string
from normalizer.position_converter import position_convert

from .commons import code_in_errors, get_errors_codes, patch_retriever


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
    assert get_errors_codes(p_c["errors"]) == {"ENOINPUTS", "ENOINPUTSOTHER"}


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
    assert get_errors_codes(p_c["errors"]) == {"ENOINPUTS", "ENOINPUTSOTHER"}


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
    assert p_c["infos"][0]["code"] == "ITOSELECTOR"


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
        include_overlapping=False,
    )
    assert model_to_string(p_c["converted_model"]) == "NG_012337.1(NM_003002.2):c.-5035"


def test_outside_cds_description():
    p_c = position_convert(
        description="NG_012337.1:27",
        to_selector_id="NM_003002.2",
        include_overlapping=False,
    )
    assert model_to_string(p_c["converted_model"]) == "NG_012337.1(NM_003002.2):c.-5035"


def test_outside_cds_description_model():
    p_c = position_convert(
        description_model=parse_description_to_model("NG_012337.1:27"),
        to_selector_id="NM_003002.2",
        include_overlapping=False,
    )
    assert model_to_string(p_c["converted_model"]) == "NG_012337.1(NM_003002.2):c.-5035"
