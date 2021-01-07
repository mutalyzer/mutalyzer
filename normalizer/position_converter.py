from mutalyzer_hgvs_parser import parse_description_to_model

from .converter import to_internal_coordinates
from .converter.to_hgvs_coordinates import to_hgvs_locations
from .description_model import location_to_description
from .position_check import check_locations
from .reference import get_reference_model
from .description import Description


class PositionConvert(object):
    def __init__(
        self,
        description="",
        reference_id="",
        position="",
        from_selector_id="",
        from_coordinate_system="",
        description_model={},
        to_selector_id="",
        to_coordinate_system="",
        include_overlapping=False,
    ):
        self.errors = []
        if description:
            self.description = Description(description=description)
            self.description_model = self.description_to_model()
        elif reference_id or position or from_selector_id or from_coordinate_system:
            self.description = Description(description_model=self._get_model_from_segmented_input(
                reference_id, position, from_selector_id, from_coordinate_system
            ))
        elif description_model:
            self.description = Description(description_model=description_model)

        else:
            self.errors.append({"code": "EINIT", "details": "Initialization error."})

        self.to_selector_id = to_selector_id
        self.to_coordinate_system = to_coordinate_system
        self.include_overlapping = include_overlapping

        self.internal_model = {}
        self.converted_model = {}

        self.description.to_internal_coordinate_model()

        self.internal_model = self.description.internal_coordinates_model

        if self.internal_model and self.to_selector_id:
            # check_locations(self.description_model, self.internal_model)
            self.converted_model = self.get_converted_model()

        self.output = self.get_output()

    def description_to_model(self):
        try:
            model = parse_description_to_model(self.description)
        except Exception as e:
            model = {
                "errors": [
                    {
                        "details": "Some error occured during description parsing.",
                        "raw_message": e,
                    }
                ]
            }
        return model

    def get_internal_model(self):
        return to_internal_coordinates.to_internal_coordinates(self.description_model, self.reference_model)

    def get_converted_model(self):
        return to_hgvs_locations(
            self.internal_model, self.description.references,
            self.to_coordinate_system, self.to_selector_id
        )

    def _get_model_from_segmented_input(
        self,
        reference_id="",
        position="",
        from_selector_id="",
        from_coordinate_system="",
    ):
        if not (reference_id and position):
            self.errors.append({"code": "ENOINPUTS"})
            return {}
        description_model = {"reference": {"id": reference_id}}

        if from_selector_id:
            description_model["reference"]["selector"] = {"id": from_selector_id}

        if from_coordinate_system:
            description_model["coordinate_system"] = from_coordinate_system

        if not isinstance(position, str):
            self.errors.append(
                {"code": "EPOSITIONINVALID", "details": "Position must be string"}
            )
            return description_model
        location_model = parse_description_to_model(position, start_rule="location")
        if location_model.get("errors"):
            self.errors.append(
                {"code": "EPOSITIONSYNTAX", "details": location_model["errors"][0]}
            )
            return description_model

        description_model["variants"] = [{"location":  location_model}]

        return description_model

    def add_overlapping(self):
        pass

    def get_output(self):
        output = {}
        if self.description.input_model:
            output["input_model"] = self.description.input_model
        if self.internal_model:
            output["internal_model"] = self.internal_model
        if self.converted_model:
            output["converted_model"] = self.converted_model
        if self.errors:
            output["errors"] = self.errors
        return output


def position_convert(
    description,
    reference_id,
    position,
    from_selector_id="",
    from_coordinate_system="",
    to_selector_id="",
    to_coordinate_system="",
    include_overlapping=False,
):
    p_c = PositionConvert(
        description=description,
        reference_id=reference_id,
        from_selector_id=from_selector_id,
        from_coordinate_system=from_coordinate_system,
        position=position,
        to_selector_id=to_selector_id,
        to_coordinate_system=to_coordinate_system,
        include_overlapping=include_overlapping,
    )
    return p_c.output
