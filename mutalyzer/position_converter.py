from mutalyzer_hgvs_parser import to_model
from mutalyzer_hgvs_parser.exceptions import UnexpectedCharacter, UnexpectedEnd

import mutalyzer.errors as errors
import mutalyzer.infos as infos

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
        self.overlap = []

        self._check_input_requirements()

        self._get_description()

        self._check_to_parameters()

        self._convert()

        self._add_overlapping()

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
            self.errors.append(errors.no_inputs())

        if not (
            self.from_selector_id
            or self.from_coordinate_system
            or self.to_selector_id
            or self.to_coordinate_system
        ):
            self.errors.append(errors.no_inputs_other())

        if self.position and not isinstance(self.position, str):
            self.errors.append(errors.position_invalid())

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
            self.description.retrieve_references()
            self.description.pre_conversion_checks()
            self.description.to_internal_indexing_model()
            self.errors += self.description.errors
            self.infos += self.description.infos

    @check_errors
    def _get_model_from_segmented_input(self):
        description_model = {"reference": {"id": self.reference_id}}

        if self.from_selector_id:
            description_model["reference"]["selector"] = {"id": self.from_selector_id}
        if self.from_coordinate_system:
            description_model["coordinate_system"] = self.from_coordinate_system

        try:
            location_model = to_model(self.position, start_rule="location")
            description_model["variants"] = [{"location": location_model}]
        except UnexpectedCharacter as e:
            self.errors.append(errors.position_syntax("Unexpected character.", e))
        except UnexpectedEnd as e:
            self.errors.append(errors.position_syntax("Unexpected end of input.", e))

        self.input_model = description_model

    @check_errors
    def _check_to_parameters(self):
        if not (self.to_selector_id or self.to_coordinate_system):
            self.to_coordinate_system = get_coordinate_system_from_reference(
                self.description.references["reference"]
            )
            if self.to_coordinate_system:
                self.infos.append(
                    infos.to_coordinate_system_from_reference(self.to_coordinate_system)
                )
            else:
                self.errors.append(
                    errors.no_to_selector(self.reference_id, self.to_selector_id)
                )

        elif self.to_selector_id:
            if not is_selector_in_reference(
                self.to_selector_id, self.description.references["reference"]
            ):
                # TODO: update the error.
                self.errors.append(
                    errors.no_to_selector(self.reference_id, self.to_selector_id)
                )

        if (self.to_selector_id == self.from_selector_id) and (
            self.to_coordinate_system == self.from_coordinate_system
        ):
            self.infos.append(infos.from_to_selector_equal())

    @check_errors
    def _convert(self):
        self.converted_model = to_hgvs_locations(
            model=self.description.internal_indexing_model,
            references=self.description.references,
            to_coordinate_system=self.to_coordinate_system,
            to_selector_id=self.to_selector_id,
            degenerate=True,
        )

    @check_errors
    def _add_overlapping(self):
        if self.include_overlapping:
            self.description.construct_equivalent(
                self.description.internal_indexing_model, False
            )
            overlap = {}
            for k in self.description.equivalent:
                if k != "g":
                    for equivalent in self.description.equivalent[k]:
                        if equivalent["reference"]["selector"]["id"] not in [
                            self.to_selector_id,
                            self.from_selector_id,
                        ]:
                            if overlap.get(k) is None:
                                overlap[k] = []
                            overlap[k].append(equivalent)
            self.overlap = overlap

    def get_output(self):
        output = {}
        if self.input_model:
            output["input_model"] = self.input_model
        if self.internal_model:
            output["internal_model"] = self.internal_model
        if self.converted_model:
            output["converted_model"] = self.converted_model
        if self.overlap:
            output["overlap"] = self.overlap
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
