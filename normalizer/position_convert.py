from mutalyzer_crossmapper import Coding, Genomic, NonCoding
from mutalyzer_hgvs_parser import parse_description_to_model

from .converter import to_hgvs, to_internal, to_internal_coordinates
from .description import location_to_description, get_errors
from .reference import (
    get_mol_type,
    get_only_selector,
    get_selector_model,
    get_selectors_overlap,
    get_reference_model,
    coordinate_system_from_mol_type,
)
from .position_check import check_positions


class PositionConvert2(object):
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
        include_overlapping=False
    ):
        self.errors = []
        if description:
            self.description = description
            self.description_model = self.description_to_model()
        elif (reference_id or position or from_selector_id or from_coordinate_system):
            self.description_model = self.get_model_from_segmented_input(
                reference_id, position, from_selector_id, from_coordinate_system)
        elif description_model:
            self.description_model = description_model
        else:
            self.errors.append({
                "code": "EINIT",
                "details": "Initialization error."})

        self.to_selector_id = to_selector_id
        self.to_coordinate_system = to_coordinate_system
        self.include_overlapping = include_overlapping

        self.internal_model = {}
        self.converted_model = {}

        if not self.errors and not get_errors(self.description_model):
            self.internal_model = self.get_internal_model()
        if self.internal_model:
            check_positions(self.description_model, self.internal_model)
            if not get_errors(self.description_model):
                self.converted_model = self.get_converted_model()

        self.output = self.get_output()

    def description_to_model(self):
        try:
            model = parse_description_to_model(self.description)
        except Exception as e:
            model = {"errors": [{
                "details": "Some error occured during description parsing.",
                "raw_message": e
            }]}
        return model

    def get_internal_model(self):
        return to_internal_coordinates.to_internal_coordinates(
            self.description_model)

    def get_converted_model(self):
        return to_internal_coordinates.to_hgvs(
            self.internal_model,
            self.to_coordinate_system,
            self.to_selector_id)

    def get_model_from_segmented_input(
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
            description_model['coordinate_system'] = from_coordinate_system

        if not isinstance(position, str):
            self.errors.append(
                {"code": "EPOSITIONINVALID", "details": "Position must be string"}
            )
            return description_model
        location_model = parse_description_to_model(position, start_rule="location")
        if location_model.get("errors"):
            self.errors.append(
                {"code": "EPOSITIONSYNTAX", "details": location_model["errors"][0]})
            return description_model

        description_model["location"] = location_model

        return description_model

    def check_internal(self):
        p = self.internal
        o = self.location_model
        seq_len = len(self.reference_model["sequence"]["seq"])
        greater = []
        smaller = []

        if p["type"] == "range":
            if p["start"]["type"] == "range":
                if p["start"]["start"].get("position"):
                    if p["start"]["start"]["position"] < 0:
                        smaller.append(location_to_description(o["start"]["start"]))
                    elif p["start"]["start"]["position"] > seq_len:
                        greater.append(location_to_description(o["start"]["start"]))
                if p["start"]["end"].get("position"):
                    if p["start"]["end"]["position"] < 0:
                        smaller.append(location_to_description(o["start"]["end"]))
                    if p["start"]["end"]["position"] > seq_len:
                        greater.append(location_to_description(o["start"]["end"]))
            elif p["start"].get("position"):
                if p["start"]["position"] < 0:
                    smaller.append(location_to_description(o["start"]))
                if p["start"]["position"] > seq_len:
                    greater.append(
                        location_to_description(o["start"]))

            if p["end"]["type"] == "range":
                if p["end"]["start"].get("position"):
                    if p["end"]["start"]["position"] < 0:
                        smaller.append(location_to_description(o["end"]["start"]))
                    if p["end"]["start"]["position"] > seq_len:
                        greater.append(location_to_description(o["end"]["start"]))
                if p["end"]["end"].get("position"):
                    if p["end"]["end"]["position"] < 0:
                        smaller.append(location_to_description(o["end"]["end"]))
                    if p["end"]["end"]["position"] > seq_len:
                        greater.append(location_to_description(o["end"]["end"]))
            elif p["end"].get("position"):
                if p["end"]["position"] <0:
                    smaller.append(location_to_description(o["end"]))
                if p["end"]["position"] > seq_len:
                    greater.append(location_to_description(o["end"]))
        elif p["type"] == "point" and p.get("position"):
            if p["position"] < 0:
                smaller.append(location_to_description(o))

            if p["position"] > seq_len:
                greater.append(location_to_description(o["position"]))
        smaller = ", ".join(smaller)
        greater = ", ".join(greater)
        if smaller:
            self.errors.append(
                {"code": "EOUTOFBOUNDARYS", "details": smaller})
        if greater:
            self.errors.append(
                {"code": "EOUTOFBOUNDARYSG", "details": greater})

    def identify_inconsistencies(self):
        if (self.to_selector_id == self.from_selector_id) and (
            self.to_coordinate_system == self.from_coordinate_system
        ):
            self.errors.append(
                {
                    "code": "EFROMTOSELECTORSEQUAL",
                    "details": "Both from and to coordinate systems are "
                    "the same, no conversion can be implemented.",
                }
            )

    def add_overlapping(self):
        if self.include_overlapping and self.internal:
            overlapping = []
            for selector in get_selectors_overlap(
                self.internal["position"], self.reference_model["model"]
            ):
                coordinate_system = coordinate_system_from_mol_type(selector)
                if coordinate_system:
                    hgvs = to_hgvs.point_to_hgvs(self.internal, **crossmap)
                if selector["id"] != self.to_selector_id:
                    overlapping.append(
                        {
                            "selector_id": selector["id"],
                            "coordinate_system": selector["coordinate_system"],
                            "position": location_to_description(hgvs),
                        }
                    )
            self.overlapping = overlapping

    def get_output(self):
        output = {}
        if self.description_model:
            output['input_model'] = self.description_model
        if self.internal_model:
            output['internal_model'] = self.internal_model
        if self.converted_model:
            output['converted_model'] = self.converted_model
        if self.errors:
            output['errors'] = self.errors
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
    p_c = PositionConvert2(
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


