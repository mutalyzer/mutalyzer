from ..reference import get_reference_model, get_mol_type, get_selector_model
from ..reference import Reference, coordinate_system_from_selector
from ..description import get_errors, get_coordinate_system, get_selector_id
from mutalyzer_crossmapper import Genomic, NonCoding, Coding


def get_point_value(point):
    value = point["position"]
    if point.get("outside_cds") and point["outside_cds"] == "upstream":
        value *= -1
    return value


def point_to_x_coding(point):
    position = point["position"]
    if point.get("outside_cds"):
        if point["outside_cds"] == "upstream":
            section = -1
            position = -1 * position
        elif point["outside_cds"] == "downstream":
            section = 1
    else:
        section = 0
    if point.get("offset"):
        offset = point["offset"]["value"]
    else:
        offset = 0
    return position, offset, section, 0


def create_exact_point_model(point):
    return {"type": "point", "position": point}


def point_to_coding(point, crossmap_function, point_function):
    if point.get("uncertain"):
        return {"type": "point", "uncertain": True}
    else:
        return create_exact_point_model(crossmap_function(point_function(point)))


def point_to_internal(point, crossmap):
    new_point = point_to_coding(point, **crossmap)
    return new_point


def range_to_internal_coordinate(location, crossmap):
    new_range = {
        "start": point_to_coding(location["start"], **crossmap),
        "end": point_to_coding(location["end"], **crossmap),
        "type": "range",
    }

    if location.get("uncertain"):
        new_range["uncertain"] = location["uncertain"]

    return new_range


def location_to_internal_coordinate(location, crossmap):
    if location["type"] == "range":
        if location["start"]["type"] == "range":
            new_location = {
                "start": range_to_internal_coordinate(location["start"], crossmap)
            }
        else:
            new_location = {"start": point_to_coding(location["start"], **crossmap)}
        if location["end"]["type"] == "range":
            new_location["end"] = range_to_internal_coordinate(
                location["end"], crossmap
            )
        else:
            new_location["end"] = point_to_coding(location["end"], **crossmap)
        new_location["type"] = "range"
        if location.get("uncertain"):
            new_location["uncertain"] = "uncertain"
    elif location["type"] == "point":
        new_location = point_to_internal(location, crossmap)

    return new_location


def variants_to_internal_coordinate(variants, crossmap):
    new_variants = []
    for variant in variants:
        import copy
        new_variant = copy.deepcopy(variant)
        if variant.get('location'):
            new_variant['location'] = location_to_internal_coordinate(
                variant["location"], crossmap)
        if variant.get('deleted'):
            for i, deleted in enumerate(variant['deleted']):
                if deleted['location']:
                    new_variant["deleted"][i]['location'] = location_to_internal_coordinate(
                deleted['location'], crossmap)
        if variant.get('inserted'):
            for i, inserted in enumerate(variant['inserted']):
                if inserted['source'] in ['reference', 'description']:
                    if inserted['location']:
                        new_variant["inserted"][i]['location'] = location_to_internal_coordinate(
                            inserted['location'], crossmap)
                else:
                    new_variant["inserted"][i] = to_internal_coordinates(inserted)
        new_variants.append(new_variant)
    return new_variants


def crossmap_setup(coordinate_system, selector_model=None):
    if coordinate_system == "g":
        crossmap = Genomic()
        return {
            "crossmap_function": crossmap.genomic_to_coordinate,
            "point_function": get_point_value,
        }
    elif coordinate_system == "c":
        crossmap = Coding(
            selector_model["exon"],
            selector_model["cds"][0],
            selector_model["inverted"],
        )
        return {
            "crossmap_function": crossmap.coding_to_coordinate,
            "point_function": point_to_x_coding,
        }
    elif coordinate_system == "n":
        crossmap = NonCoding(selector_model["exon"], selector_model["inverted"])
        return {
            "crossmap_function": crossmap.noncoding_to_coordinate,
            "point_function": point_to_x_coding,
        }


def add_msg(dictionary, message_type, message):
    if dictionary.get(message_type) is None:
        dictionary[message_type] = []
    dictionary[message_type].append(message)


def get_crossmapper_inputs_from_selector_id(description_model, reference):
    selector_id = get_selector_id(description_model)
    selector_model = reference.get_selector_model(selector_id)
    if selector_model:
        selector_coordinate_system = coordinate_system_from_selector(
            selector_model)
        if selector_coordinate_system:
            coordinate_system = selector_coordinate_system
            description_model['coordinate_system'] = coordinate_system
            add_msg(description_model, 'info', {
                "code": "ICOORDINATESYSTEM",
                "details": "Coordinate system identified as {} "
                           "from the reference molecule type.".format(
                    selector_coordinate_system)})
            return coordinate_system, selector_model
        else:
            add_msg(description_model, 'errors', {"code": "ENOCOORDINATESYSTEM",
                  "details": "No coordinate system provided and it "
                             "cannot be identified from the {} "
                             "selector.".format(selector_id)})
            return None, None
    else:
        add_msg(description_model, 'errors', {"code": "ENOCOORDINATESYSTEM",
              "details": "No coordinate system provided and it "
                         "cannot be identified from the {} "
                         "selector.".format(selector_id)})
        add_msg(description_model['reference']['selector'], 'errors',
                {"code": "ENOSELECTORID",
              "details": "No {} selector found in reference {}.".format(
                  selector_id, reference.id)})
        return None, None


def get_crossmapper_inputs(description_model, reference):
    coordinate_system = get_coordinate_system(description_model)

    if not coordinate_system:
        # Try to identify a coordinate system.
        # 1. From the selector id.
        if get_selector_id(description_model):
            return get_crossmapper_inputs_from_selector_id(
                description_model, reference)
        # 2. From the reference.
        reference_coordinate_system = reference.get_default_coordinate_system()
        if reference_coordinate_system:
            coordinate_system = reference_coordinate_system
            description_model['coordinate_system'] = coordinate_system
            add_msg(description_model, 'info', {
                "code": "ICOORDINATESYSTEM",
                "details": "Coordinate system identified as {} "
                           "from the reference molecule type.".format(
                    reference_coordinate_system),
                })
            # Do not return yet, we need to see if we also have a model.
        else:
            add_msg(description_model['reference'], 'errors', {"code": "ENOCOORDINATESYSTEM",
                  "details": "No coordinate system provided and it "
                             "cannot be identified from the reference."})
            return None, None

    # We have a coordinate system, either provided or determined from the
    # reference. We need to check if it is valid and if we also have a selector
    # model.
    selector_id = get_selector_id(description_model)
    if selector_id:
        selector_model = reference.get_selector_model(selector_id)
        if selector_model:
            selector_coordinate_system = coordinate_system_from_selector(
                selector_model)
            if coordinate_system == selector_coordinate_system:
                return coordinate_system, selector_model
            else:
                add_msg(description_model, 'errors', {"code": "ECOORDINATESYSTEM",
                      "details": "Coordinate system {} does not match with {} "
                                 "selector {} coordinate system. ".format(
                          coordinate_system, selector_id, selector_coordinate_system)})
                return None, None
        else:
            add_msg(description_model['reference']['selector'], 'errors',
                    {"code": "ENOSELECTORID",
                  "details": "No {} selector found in reference {}.".format(
                      selector_id, reference.id)})
            return None, None

    if coordinate_system in ['g']:
        return coordinate_system, None

    # For some coordinates we still need a selector model.
    if coordinate_system in ['c', 'n']:
        only_selector = reference.get_only_selector()
        if only_selector:
            add_msg(description_model, 'info', {"code": "IONLYSELECTOR",
                "details": "Selector {} identified as the only one in the reference.".format(
                    only_selector['id']),
                })
            return coordinate_system, only_selector
        else:
            add_msg(description_model['reference'], 'errors', {"code": "ESELECTORREQUIRED",
                  "details": "Selector not provided, but required for the {} "
                             "coordinate system, or it should be in reference "
                             "and its not.".format(coordinate_system)})
            return None, None


def to_internal_coordinates(description_model):
    print(description_model)
    if description_model.get('reference') and description_model['reference'].get('id'):
        reference = Reference(description_model['reference']['id'])
        if not reference.model:
            add_msg(description_model,
                    'errors', {"code": "ERETR", "details": "Reference not retrieved"})
            # Without a reference we cannot proceed further.
            return
    else:
        if description_model.get('reference') and description_model['reference'].get('id'):
        add_msg(description_model,
                'errors',
                {"code": "ERETR", "details": "Reference not retrieved"})
        # Without a reference we cannot proceed further.
        return

        return description_model['errors'] = [{"code": "ERETR"}]

    # We have a reference model, we can proceed further
    coordinate_system, selector_model = get_crossmapper_inputs(description_model, reference)
    if coordinate_system:
        if get_errors(description_model):
            print("Something strange happened.")
        crossmap = crossmap_setup(coordinate_system, selector_model)
        internal_variants = variants_to_internal_coordinate(description_model['variants'], crossmap)
        print('internal variants:')
        print(internal_variants)
