from mutalyzer_crossmapper import Coding, Genomic, NonCoding
from mutalyzer_hgvs_parser import parse_description_to_model

from .converter import to_hgvs, to_internal
from .description import location_to_description
from .normalizer import get_reference_model
from .reference import (
    get_mol_type,
    get_only_selector,
    get_selector_model,
    get_selectors_overlap,
)


def get_coordinate_system(selector_model):
    if selector_model["type"] in ["mRNA"]:
        return "c"
    elif selector_model["type"] in ["ncRNA"]:
        return "n"
    else:
        return ""


def crossmap_to_x_setup(coordinate_system, selector_model=None):
    if coordinate_system == "g":
        crossmap = Genomic()
        return {
            "crossmap_function": crossmap.genomic_to_coordinate,
            "point_function": to_internal.get_point_value,
        }
    elif coordinate_system == "c":
        crossmap = Coding(
            selector_model["exon"],
            selector_model["cds"][0],
            selector_model["inverted"],
        )
        return {
            "crossmap_function": crossmap.coding_to_coordinate,
            "point_function": to_internal.point_to_x_coding,
        }
    elif coordinate_system == "n":
        crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
        return {
            "crossmap_function": crossmap.noncoding_to_coordinate,
            "point_function": to_internal.point_to_x_coding,
        }


def crossmap_to_hgvs_setup(coordinate_system, selector_model):
    if coordinate_system == "g":
        crossmap = Genomic()
        return {
            "crossmap_function": crossmap.coordinate_to_genomic,
            "point_function": to_hgvs.genomic_to_point,
        }
    elif coordinate_system == "c":
        crossmap = Coding(
            selector_model["exon"],
            selector_model["cds"][0],
            selector_model["inverted"],
        )
        return {
            "crossmap_function": crossmap.coordinate_to_coding,
            "point_function": to_hgvs.coding_to_point,
            "degenerate": True,
        }
    elif coordinate_system == "n":
        crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
        return {
            "crossmap_function": crossmap.coordinate_to_noncoding,
            "point_function": to_hgvs.noncoding_to_point,
        }


class PositionConvert(object):
    def __init__(
        self,
        reference_id,
        position,
        from_selector_id="",
        from_coordinate_system="",
        to_selector_id="",
        to_coordinate_system="",
        include_overlapping=False,
    ):
        self.reference_id = reference_id
        self.position = position
        self.from_selector_id = from_selector_id
        self.from_coordinate_system = from_coordinate_system
        self.to_selector_id = to_selector_id
        self.to_coordinate_system = to_coordinate_system
        self.include_overlapping = include_overlapping

        self.reference_model = {}
        self.from_selector_model = {}
        self.to_selector_model = {}
        self.mol_type = ""
        self.reference_coordinate_system = ""
        self.location_model = {}

        self.internal = None
        self.converted = None
        self.overlapping = []

        self.errors = []
        self.infos = []

        self.output = {}

        self.process_inputs()

        self.convert()

        self.add_overlapping()

        self.construct_output()

    def convert(self):
        if self.errors:
            return
        crossmap = crossmap_to_x_setup(
            self.from_coordinate_system, self.from_selector_model
        )
        self.internal = to_internal.point_to_coding(self.location_model, **crossmap)

        if self.internal["position"] < 0:
            self.errors.append({"code": "EOUTOFBOUNDARY"})
        if self.internal["position"] > len(self.reference_model["sequence"]["seq"]):
            self.errors.append({"code": "EOUTOFBOUNDARY"})

        if self.errors:
            return
        crossmap = crossmap_to_hgvs_setup(
            self.to_coordinate_system, self.to_selector_model
        )
        self.converted = location_to_description(
            to_hgvs.point_to_hgvs(self.internal, **crossmap)
        )

    def process_inputs(self):
        """
        Calls the appropriate input processing functions to generate a summary
        of the inputs.
        """
        if not (self.reference_id and self.position):
            self.errors.append({"code": "ENOINPUTS"})
            return

        if not (
            self.from_selector_id
            or self.from_coordinate_system
            or self.to_selector_id
            or self.to_coordinate_system
        ):
            self.errors.append({"code": "ENOINPUTSOTHER"})
            return

        self.process_reference_id()
        self.process_from_selector_id()
        self.process_from_coordinate_system()
        self.process_to_selector_id()
        self.process_to_coordinate_system()
        self.identify_inconsistencies()
        self.process_position()

    def process_reference_id(self):
        if self.reference_id:
            self.reference_model = get_reference_model(self.reference_id)
            if self.reference_model:
                self.mol_type = get_mol_type(self.reference_model)
                if self.mol_type in ["genomic DNA", "dna"]:
                    self.reference_coordinate_system = "g"
                else:
                    self.errors.append(
                        {
                            "code": "EUNSUPPORTEDREF",
                            "details": "Reference {} mol_type not supported.".format(
                                self.mol_type
                            ),
                        }
                    )
            else:
                self.errors.append({"code": "ERETR"})

    def process_from_selector_id(self):
        if self.from_selector_id and self.reference_model:
            self.from_selector_model = get_selector_model(
                self.reference_model["model"], self.from_selector_id
            )
            if not self.from_selector_model:
                self.errors.append({"code": "ENOFROMSELECTOR"})

    def process_from_coordinate_system(self):
        if self.from_coordinate_system == "g":
            if self.from_selector_id:
                self.errors.append(
                    {
                        "code": "EINCONSISTENINPUTS",
                        "details": "Genomic coordinate system with selector id "
                        "provided, which one to choose?",
                    }
                )
            elif self.mol_type not in ["genomic DNA", "dna"]:
                self.errors.append(
                    {
                        "code": "EFROMSELECTORCS",
                        "details": "Coordinate system does not match the reference.",
                    }
                )
        elif self.from_coordinate_system == "c":
            if self.from_selector_id and self.from_selector_model:
                if self.from_selector_model["type"] not in ["mRNA"]:
                    self.errors.append(
                        {
                            "code": "EFROMSELECTORCS",
                            "details": "Coordinate system does not match the selector.",
                        }
                    )
            else:
                self.try_only_selector('from')
        elif self.from_coordinate_system == "n":
            if self.from_selector_id and self.from_selector_model:
                if self.from_selector_model["type"] not in ["ncRNA"]:
                    self.errors.append(
                        {
                            "code": "EFROMSELECTORCS",
                            "details": "Coordinate system does not match the selector.",
                        }
                    )
            else:
                self.try_only_selector('from')
        elif self.from_coordinate_system == "Selector" or (
            self.from_coordinate_system == "" and self.from_selector_id
        ):
            if not self.from_selector_id:
                self.try_only_selector('from')
                # self.errors.append(
                #     {
                #         "code": "EFROMSELECTORCS",
                #         "details": "Selector id must be provided in order to "
                #         "identify its coordinate system.",
                #     }
                # )
            elif self.from_selector_model:
                self.from_coordinate_system = get_coordinate_system(
                    self.from_selector_model
                )
                if self.from_coordinate_system:
                    self.infos.append(
                        {
                            "code": "IFROMSELECTOR",
                            "details": "From coordinate system identified as {} "
                            "from the selector molecule type.".format(
                                self.from_coordinate_system
                            ),
                        }
                    )
        elif self.from_coordinate_system in ["Reference", ""] or not (
            self.from_coordinate_system
        ):
            if self.mol_type in ["genomic DNA", "dna"]:
                self.from_coordinate_system = "g"
                self.infos.append(
                    {
                        "code": "IFROMSELECTOR",
                        "details": "from_coordinate_system identified as g "
                        "from the reference molecule type.",
                    }
                )
            # else: we should not reach, since the reference is checked first.

    def process_to_selector_id(self):
        if self.to_selector_id and self.reference_model:
            self.to_selector_model = get_selector_model(
                self.reference_model["model"], self.to_selector_id
            )
            if not self.to_selector_model:
                self.errors.append({"code": "ENOTOSELECTOR"})

    def process_to_coordinate_system(self):
        if self.to_coordinate_system == "g":
            if self.to_selector_id:
                self.errors.append(
                    {
                        "code": "EINCONSISTENINPUTS",
                        "details": "Genomic coordinate system with selector id "
                        "provided, which one to choose?",
                    }
                )
            elif self.mol_type not in ["genomic DNA", "dna"]:
                self.errors.append(
                    {
                        "code": "ETOSELECTORCS",
                        "details": "Coordinate system does not match the reference.",
                    }
                )
        elif self.to_coordinate_system == "c":
            if self.to_selector_id and self.to_selector_model:
                if self.to_selector_model["type"] not in ["mRNA"]:
                    self.errors.append(
                        {
                            "code": "ETOSELECTORCS",
                            "details": "Coordinate system does not match the selector.",
                        }
                    )
            else:
                self.try_only_selector('to')
        elif self.to_coordinate_system == "n":
            if self.to_selector_id and self.to_selector_model:
                if self.to_selector_model["type"] not in ["ncRNA"]:
                    self.errors.append(
                        {
                            "code": "ETOSELECTORCS",
                            "details": "Coordinate system does not match the selector.",
                        }
                    )
            else:
                self.try_only_selector('to')
        elif self.to_coordinate_system == "Selector" or (
            self.to_coordinate_system == "" and self.to_selector_id
        ):
            if not self.to_selector_id:
                self.try_only_selector('to')
                # self.errors.append(
                #     {
                #         "code": "ETOSELECTORCS",
                #         "details": "Selector id must be provided in order to "
                #         "identify its coordinate system.",
                #     }
                # )
            elif self.to_selector_model:
                self.to_coordinate_system = get_coordinate_system(
                    self.to_selector_model
                )
                if self.to_coordinate_system:
                    self.infos.append(
                        {
                            "code": "ITOSELECTOR",
                            "details": "To coordinate system identified as {} "
                            "from the selector molecule type.".format(
                                self.to_coordinate_system
                            ),
                        }
                    )
        elif (
            self.to_coordinate_system in ["Reference", ""]
            or not self.to_coordinate_system
        ):
            if self.mol_type in ["genomic DNA", "dna"]:
                self.to_coordinate_system = "g"
                self.infos.append(
                    {
                        "code": "ITOSELECTOR",
                        "details": "to_coordinate_system identified as g "
                        "from the reference molecule type.",
                    }
                )
            # else: we should not reach, since the reference is checked first.

    def try_only_selector(self, source):
        if source == 'from':
            if self.from_coordinate_system == 'Selector':
                only_selector = get_only_selector(self.reference_model["model"])
            else:
                only_selector = get_only_selector(
                    self.reference_model["model"], self.from_coordinate_system)
        elif source == 'to':
            if self.to_coordinate_system == 'Selector':
                only_selector = get_only_selector(self.reference_model["model"])
            else:
                only_selector = get_only_selector(
                    self.reference_model["model"], self.to_coordinate_system)

        if only_selector:
            if source == 'from':
                self.from_selector_id = only_selector["id"]
                self.from_selector_model = only_selector
                self.infos.append({"code": "IFROMONLYSELECTOR"})
                if self.from_coordinate_system == 'Selector':
                    self.from_coordinate_system = get_coordinate_system(only_selector)
                    self.infos.append({"code": "IFROMONLYCOORDINATESYSTEM"})
            elif source == "to":
                self.to_selector_id = only_selector["id"]
                self.to_selector_model = only_selector
                self.infos.append({"code": "ITOONLYSELECTOR"})
                if self.to_coordinate_system == 'Selector':
                    self.to_coordinate_system = get_coordinate_system(only_selector)
                    self.infos.append({"code": "ITOONLYCOORDINATESYSTEM"})
        else:
            if source == 'from':
                self.errors.append({"code": "ENOFROMSELECTOR"})
            elif source == 'to':
                self.errors.append({"code": "ENOTOSELECTOR"})

    def process_position(self):
        if not isinstance(self.position, str):
            self.errors.append(
                {"code": "EPOSITIONINVALID", "details": "Position must be string"}
            )
            return
        self.location_model = parse_description_to_model(
            self.position, start_rule="location"
        )
        if self.location_model.get("errors"):
            self.errors.append(
                {"code": "ESYNTAX", "details": self.location_model["errors"][0]}
            )
        if self.location_model["type"] == "range":
            self.errors.append({"code": "ERANGELOCATION"})

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
            other_selectors = []
            for selector in get_selectors_overlap(
                self.internal["position"], self.reference_model["model"]
            ):
                crossmap = crossmap_to_hgvs_setup(
                    get_coordinate_system(selector), selector
                )
                hgvs = to_hgvs.point_to_hgvs(self.internal, **crossmap)
                if selector["id"] != self.to_selector_id:
                    other_selectors.append(
                        {
                            "selector_id": selector["id"],
                            "coordinate_system": selector["coordinate_system"],
                            "position": location_to_description(hgvs),
                        }
                    )
            self.overlapping = other_selectors

    def construct_output(self):
        if self.errors:
            self.output["errors"] = self.errors
        if self.infos:
            self.output["infos"] = self.infos

        if self.errors:
            return

        if self.converted and self.internal:
            input_position = {
                "coordinate_system": self.from_coordinate_system,
                "position": self.position,
            }
            if self.from_selector_id:
                input_position["selector_id"] = self.from_selector_id
            converted_position = {
                "coordinate_system": self.to_coordinate_system,
                "position": self.converted,
            }
            if self.to_selector_id:
                converted_position["selector_id"] = self.to_selector_id
            self.output["conversion"] = {
                "reference_id": self.reference_id,
                "input_position": input_position,
                "converted_position": converted_position,
            }
            if self.include_overlapping:
                self.output["conversion"]["overlapping"] = self.overlapping


def position_convert(
    reference_id,
    position,
    from_selector_id="",
    from_coordinate_system="",
    to_selector_id="",
    to_coordinate_system="",
    include_overlapping=False,
):
    p_c = PositionConvert(
        reference_id=reference_id,
        from_selector_id=from_selector_id,
        from_coordinate_system=from_coordinate_system,
        position=position,
        to_selector_id=to_selector_id,
        to_coordinate_system=to_coordinate_system,
        include_overlapping=include_overlapping,
    )
    return p_c.output
