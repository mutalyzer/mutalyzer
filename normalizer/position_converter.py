from mutalyzer_hgvs_parser import parse_description_to_model

from .converter.to_hgvs_coordinates import to_hgvs_locations
from .description import Description
from .reference import get_coordinate_system_from_reference, is_selector_in_reference
from .util import check_errors


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
        print("-----")
        self.input_description = description
        self.reference_id = reference_id
        self.position = position
        self.from_selector_id = from_selector_id
        self.from_coordinate_system = from_coordinate_system
        self.input_model = description_model
        self.to_selector_id = to_selector_id
        self.to_coordinate_system = to_coordinate_system
        self.include_overlapping = include_overlapping

        self.stop_on_errors = False

        self.errors = []
        self.infos = []

        self.description = None

        self.internal_model = {}
        self.converted_model = {}

        self._check_input_requirements()

        self._get_description()

        self._check_to_parameters()

        self._convert()

        self.output = self.get_output()

    def _check_input_requirements(self):
        """
        Checks if the provided inputs can be used for conversion.
        """
        if not (
            self.input_description
            or self.input_model
            or (self.reference_id and self.position)
        ):
            self.errors.append({"code": "ENOINPUTS"})

        if not (
            self.from_selector_id
            or self.from_coordinate_system
            or self.to_selector_id
            or self.to_coordinate_system
        ):
            self.errors.append({"code": "ENOINPUTSOTHER"})

    @check_errors
    def _get_description(self):
        if self.input_description:
            self.description = Description(description=self.input_description)
        elif self.input_model:
            self.description = Description(description_model=self.input_model)
        elif (
            self.reference_id
            or self.position
            or self.from_selector_id
            or self.from_coordinate_system
        ):
            self._get_model_from_segmented_input()
            # TODO: Check the schema.
            self.description = Description(description_model=self.input_model)

        if self.description:
            self.description.to_internal_indexing_model()
            self.errors += self.description.errors
            self.infos += self.description.infos

    @check_errors
    def _get_model_from_segmented_input(self):
        description_model = {"reference": {"id": self.reference_id}}

        if self.from_selector_id:
            description_model["reference"]["selector"] = {"id": self.from_selector_id}

        if (
            self.from_coordinate_system
            and self.from_coordinate_system.strip().lower() not in ["selector"]
        ):
            description_model["coordinate_system"] = self.from_coordinate_system

        if not isinstance(self.position, str):
            self.errors.append(
                {"code": "EPOSITIONINVALID", "details": "Position must be string"}
            )
            return
        location_model = parse_description_to_model(
            self.position, start_rule="location"
        )
        if location_model.get("errors"):
            self.errors.append(
                {"code": "EPOSITIONSYNTAX", "details": location_model["errors"][0]}
            )
            return

        description_model["variants"] = [{"location": location_model}]

        self.input_model = description_model

    @check_errors
    def _check_to_parameters(self):
        if not (self.to_selector_id or self.to_coordinate_system):
            self.to_coordinate_system = get_coordinate_system_from_reference(
                self.description["references"]["reference"]
            )
        elif self.to_selector_id:
            if not is_selector_in_reference(
                self.to_selector_id, self.description.references["reference"]
            ):
                # TODO: update the error.
                self.errors.append({"code": "ENOTOSELECTOR"})

    @check_errors
    def _convert(self):
        self.converted_model = to_hgvs_locations(
            model=self.description.internal_indexing_model,
            references=self.description.references,
            to_coordinate_system=self.to_coordinate_system,
            to_selector_id=self.to_selector_id,
            degenerate=True,
        )

    def add_overlapping(self):
        pass

    def get_output(self):
        output = {}
        if self.input_model:
            output["input_model"] = self.input_model
        if self.internal_model:
            output["internal_model"] = self.internal_model
        if self.converted_model:
            output["converted_model"] = self.converted_model
        if self.errors:
            output["errors"] = self.errors
        if self.infos:
            output["infos"] = self.infos
        return output


def position_convert(
    description="",
    reference_id="",
    position="",
    from_selector_id="",
    from_coordinate_system="",
    description_model=None,
    to_selector_id="",
    to_coordinate_system="",
    include_overlapping=False,
):
    p_c = PositionConvert(
        description=description,
        reference_id=reference_id,
        position=position,
        from_selector_id=from_selector_id,
        from_coordinate_system=from_coordinate_system,
        description_model=description_model,
        to_selector_id=to_selector_id,
        to_coordinate_system=to_coordinate_system,
        include_overlapping=include_overlapping,
    )
    return p_c.output
