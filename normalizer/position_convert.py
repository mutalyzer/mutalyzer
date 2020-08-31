from mutalyzer_crossmapper import Coding, Genomic, NonCoding
from mutalyzer_hgvs_parser import parse_description_to_model

from .converter import to_hgvs, to_internal
from .description import location_to_description
from .normalizer import get_reference_model
from .reference import get_mol_type, get_selector_model, get_selectors_overlap


def crossmap_to_x_setup(selector_model, mol_type, relative_to):
    if mol_type in ["genomic DNA", "dna"]:
        if relative_to == "Reference":
            crossmap = Genomic()
            return {
                "crossmap_function": crossmap.genomic_to_coordinate,
                "point_function": to_internal.get_point_value,
            }
        elif relative_to == "Selector" and selector_model["type"] in ["mRNA"]:
            crossmap = Coding(
                selector_model["exon"],
                selector_model["cds"][0],
                selector_model["inverted"],
            )
            return {
                "crossmap_function": crossmap.coding_to_coordinate,
                "point_function": to_internal.point_to_x_coding,
            }
        elif relative_to == "Selector" and selector_model["type"] in [
            "lnc_RNA",
            "ncRNA",
        ]:
            crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
            return {
                "crossmap_function": crossmap.noncoding_to_coordinate,
                "point_function": to_internal.point_to_x_coding,
            }


def crossmap_to_hgvs_setup(selector_model, mol_type, relative_to):
    if relative_to == "Reference":
        if selector_model["type"] in ["mRNA"]:
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
        if selector_model["type"] in ["lnc_RNA", "ncRNA"]:
            crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
            return {
                "crossmap_function": crossmap.coordinate_to_noncoding,
                "point_function": to_hgvs.noncoding_to_point,
            }
    elif relative_to == "Selector":
        if mol_type in ["genomic DNA", "dna"]:
            crossmap = Genomic()
            return {
                "crossmap_function": crossmap.coordinate_to_genomic,
                "point_function": to_hgvs.genomic_to_point,
            }


def position_convert_(
    reference_id, position, from_selector_id=None, from_coordinate_system=None,
        to_selector_id=None, to_coordinate_system=None, include_overlapping=False
):
    # Retrieve the reference model. Identify its mol_type and check if is within
    # the supported ones. Set the reference coordinate system. Report any errors.
    reference_model = get_reference_model(reference_id)
    if reference_model:
        mol_type = get_mol_type(reference_model)
        if mol_type in ["genomic DNA", "dna"]:
            reference_coordinate_system = "g"
        else:
            return {"errors": [{"code": "EMOLTYPE", "details": mol_type}]}
    else:
        return {"errors": [{"code": "ERETR"}]}

    # Retrieve the location model and report any syntax errors.
    location_model = parse_description_to_model(position, start_rule="location")
    if location_model.get("errors"):
        return {"errors": [{"code": "ESYNTAX", "details": location_model["errors"][0]}]}
    if location_model["type"] == "range":
        return {"errors": [{"code": "ERANGELOCATION"}]}

    if from_selector_id is None and from_coordinate_system is None:
        from_coordinate_system = reference_coordinate_system
    # If from_selector_id supplied, try to retrieve the selector model. However,
    # maybe is better to check whether the
    if from_selector_id:
        selector_model = get_selector_model(reference_model["model"], from_selector_id)
    if selector_model:
        if selector_model["type"] in ["mRNA"]:
            selector_coordinate_system = "c"
        elif selector_model["type"] in ["lnc_RNA", "ncRNA"]:
            selector_coordinate_system = "n"
    else:
        return {"errors": [{"code": "ENOSELECTOR"}]}
    crossmap = crossmap_to_x_setup(selector_model, mol_type, relative_to)
    internal = to_internal.point_to_coding(location_model, **crossmap)

    if internal["position"] < 0:
        return {"errors": [{"code": "EOUTOFBOUNDARY"}]}
    if internal["position"] > len(reference_model["sequence"]["seq"]):
        return {"errors": [{"code": "EOUTOFBOUNDARY"}]}

    crossmap = crossmap_to_hgvs_setup(selector_model, mol_type, relative_to)
    hgvs = to_hgvs.point_to_hgvs(internal, **crossmap)
    if relative_to == "Reference":
        output = {
            "reference": {
                "id": reference_id,
                "position": position,
                "coordinate_system": reference_coordinate_system,
            },
            "selector": {
                "id": from_selector_id,
                "position": location_to_description(hgvs),
                "coordinate_system": selector_coordinate_system,
            },
        }
    elif relative_to == "Selector":
        output = {
            "reference": {
                "id": reference_id,
                "position": location_to_description(hgvs),
                "coordinate_system": reference_coordinate_system,
            },
            "selector": {
                "id": from_selector_id,
                "position": position,
                "coordinate_system": selector_coordinate_system,
            },
        }
    if include_overlapping:
        output["other_selectors"] = []
        for selector in get_selectors_overlap(
            internal["position"], reference_model["model"]
        ):
            crossmap = crossmap_to_hgvs_setup(selector, mol_type, "Reference")
            hgvs = to_hgvs.point_to_hgvs(internal, **crossmap)
            if selector["id"] != from_selector_id:
                output["other_selectors"].append(
                    {
                        "id": selector["id"],
                        "position": location_to_description(hgvs),
                        "coordinate_system": selector["coordinate_system"],
                    }
                )
    return output



class PositionConvert(object):
    def __init__(self, reference_id, position, from_selector_id=None,
                 from_coordinate_system=None, to_selector_id=None,
                 to_coordinate_system=None, include_overlapping=False):
        self.reference_id = reference_id
        self.position = position
        self.from_selector_id = from_selector_id
        self.from_coordinate_system = from_coordinate_system
        self.to_selector_id = to_selector_id
        self.to_coordinate_system = to_coordinate_system

        self.reference_model = {}
        self.from_selector_model = {}
        self.to_selector_model = {}
        self.mol_type = ''
        self.reference_coordinate_system = ''

        self.errors = []
        self.infos = []
        self.summary = {}

        self.process_inputs()
        self.update_status()

    def update_status(self):
        pass

    def ongoing(self):
        if self.errors:
            return False
        else:
            return True

    def status_check(func):
        def _decorator(self, *args, **kwargs):
            print(self.errors)

            print('self is %s' % self)
            return func(self, *args, **kwargs)

        return _decorator

    def process_inputs(self):
        """
        Calls the appropriate input processing functions to generate a summary
        of the inputs.
        """
        if not (self.reference_id and self.position):
            self.errors.append({"code": "ENOINPUTS"})

        self.process_reference_id()
        self.process_from_selector_id()
        self.process_from_coordinate_system()
        self.process_to_selector_id()

    def process_reference_id(self):
        if self.reference_id:
            self.reference_model = get_reference_model(self.reference_id)
            if self.reference_model:
                self.mol_type = get_mol_type(self.reference_model)
                if self.mol_type in ["genomic DNA", "dna"]:
                    self.reference_coordinate_system = "g"
                else:
                    self.errors.append({"code": "EMOLTYPE", "details": self.mol_type})
            else:
                self.errors.append({"code": "ERETR"})

    def process_from_selector_id(self):
        if self.from_selector_id and self.reference_model:
            self.from_selector_model = get_selector_model(
                self.reference_model["model"],
                self.from_selector_id)
            if self.from_selector_model is None:
                self.errors.append({"code": "ENOFROMSELECTOR"})

    def process_from_coordinate_system(self):
        if self.from_coordinate_system and self.from_selector_id and self.from_selector_model:
            if self.from_coordinate_system == 'c':
                if self.from_selector_model["type"] not in ["mRNA"]:
                    self.errors.append({"code": "EFROMSELECTORCS"})
            if self.from_coordinate_system == 'n':
                if self.from_selector_model["type"] not in ["ncRNA"]:
                    self.errors.append({"code": "EFROMSELECTORCS"})
        elif self.from_coordinate_system:
            if self.from_coordinate_system == 'g':
                if self.mol_type not in ["genomic DNA", "dna"]:
                    self.errors.append({"code": "EFROMSELECTORCS"})
            elif self.from_coordinate_system == 'Reference':
                if self.mol_type in ["mRNA"]:
                    self.from_coordinate_system = 'c'
                    self.infos.append(
                        {"code": "IFROMSELECTOR",
                         "details": "From coordinate system identified as c "
                                    "from the reference molecule type."})
                elif self.mol_type in ["genomic DNA", "dna"]:
                    self.from_coordinate_system = 'g'
                    self.infos.append(
                        {"code": "IFROMSELECTOR",
                         "details": "from_coordinate_system identified as g "
                                    "from the reference molecule type."})
            elif self.from_coordinate_system == 'Selector':
                if self.from_selector_id is None:
                    self.errors.append({"code": "EFROMSELECTORCS"})
                else:
                    print("yes")
                    if self.from_selector_model:
                        if self.from_selector_model["mol_type"] in ["mRNA"]:
                            self.from_coordinate_system = 'c'
                            self.infos.append(
                                {"code": "IFROMSELECTOR",
                                 "details": "From coordinate system identified as c "
                                            "from the selector molecule type."})

    def process_position(self):
        location_model = parse_description_to_model(self.position,
                                                    start_rule="location")
        if location_model.get("errors"):
            return {"errors": [
                {"code": "ESYNTAX", "details": location_model["errors"][0]}]}
        if location_model["type"] == "range":
            return {"errors": [{"code": "ERANGELOCATION"}]}

    def process_to_selector_id(self):
        if self.to_selector_id and self.reference_model:
            self.to_selector_model = get_selector_model(
                self.reference_model["model"],
                self.to_selector_id)
            if self.to_selector_model is None:
                self.errors.append({"code": "ENOTOSELECTOR"})


def position_convert(
    reference_id, position, from_selector_id=None, from_coordinate_system=None,
        to_selector_id=None, to_coordinate_system=None, include_overlapping=False
):
    p_c = PositionConvert(
        reference_id=reference_id,
        from_selector_id=from_selector_id,
        from_coordinate_system=from_coordinate_system,
        position=position,
        to_selector_id=to_selector_id,
        to_coordinate_system=to_coordinate_system,
        include_overlapping=include_overlapping)

    return {"errors": p_c.errors,
            "infos": p_c.infos}

