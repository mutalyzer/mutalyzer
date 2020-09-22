from .converter.to_internal_coordinates import get_reference
from .description import point_to_description
from .util import sort_variants, add_msg


def check_location_sequence_boundary(hgvs_location, internal_location, reference):
    if hgvs_location['type'] == 'range':
        check_location_sequence_boundary(
            hgvs_location['start'], internal_location['start'], reference)
        check_location_sequence_boundary(
            hgvs_location['end'], internal_location['end'], reference)
    if hgvs_location['type'] == 'point' and not hgvs_location.get("uncertain"):
        if (internal_location['position'] > reference.get_length()
                or internal_location['position'] < 0):
            add_msg(hgvs_location, "errors", {
                "code": "EOUTOFBOUNDARY",
                "details": "Position {} is out of sequence boundaries.".format(
                    point_to_description(hgvs_location))})
            internal_location['errors'] = [{
                "code": "EOUTOFBOUNDARY",
                "details": "Position {} is out of sequence boundaries.".format(
                    point_to_description(internal_location))
            }]


def check_points(hgvs_model, internal_model, reference):
    for k in hgvs_model.keys():
        if k in ['variants', 'deleted', 'inserted']:
            if isinstance(hgvs_model[k], list):
                for i, v in enumerate(hgvs_model[k]):
                    check_points(
                        hgvs_model[k][i], internal_model[k][i], reference)
            elif isinstance(hgvs_model[k], dict):
                check_points(
                    hgvs_model[k], internal_model[k], reference)
        elif k == 'location':
            check_location_sequence_boundary(hgvs_model[k], internal_model[k], reference)
            check_start_end(hgvs_model[k], internal_model[k])


def is_start_end_range_ok(hgvs_location, internal_location):
    if internal_location["start"].get("uncertain") or internal_location["end"].get("uncertain"):
        # ?_20, or ?_? or 20_?, we don't care.
        return False
    if internal_location["start"]["position"] > internal_location["end"]["position"]:
        # This is clear, e.g., 30_20.
        add_msg(hgvs_location, "errors", {"code": "ERANGE",
                                     "details": "Start location greater than end location."})
        return False
    return True


def check_start_end(hgvs_location, internal_location):
    if internal_location["type"] == "point":
        # There is nothing to check here.
        return

    # We have a range location, so we can start the checking.
    if internal_location["start"]["type"] == internal_location["end"]["type"] == "point":
        is_start_end_range_ok(hgvs_location, internal_location)
    elif internal_location["start"]["type"] == "range" and internal_location["end"]["type"] == "point":
        if not is_start_end_range_ok(hgvs_location["start"], internal_location["start"]):
            return
        elif internal_location["end"].get("uncertain"):
            return
        # elif location["start"]


def yield_all_locations(model):
    """
    Yields all the location models present in the input, at all their levels.
    For (10_?)_30 there will be models for (10_?)_30, (10_?), 10, ?, and 30.

    :param model: description model
    """
    for k in model.keys():
        if k in ['variants', 'deleted', 'inserted']:
            if isinstance(model[k], list):
                for sub_model in model[k]:
                    yield from yield_all_locations(sub_model)
            elif isinstance(model[k], dict):
                yield from yield_all_locations(model[k])
        elif k == 'location':
            yield model[k]
            if model[k]["type"] == "range":
                yield model[k]["start"]
                if model[k]["start"]["type"] == "range":
                    yield model[k]["start"]["start"]
                    yield model[k]["start"]["end"]
                if model[k]["end"]["type"] == "range":
                    yield model[k]["end"]["start"]
                    yield model[k]["end"]["end"]
                yield model[k]["end"]


def yield_locations_with_path(model, path=[]):
    for k in model.keys():
        if k == 'variants':
            path.append('variants')
            for i, v in enumerate(model[k]):
                path.append(i)
                yield from yield_locations_with_path(v, path)
                path.pop()
            path.pop()
        elif k == 'inserted':
            for i, v in enumerate(model[k]):
                if v.get('source') in ['description', 'reference']:
                    path.append(i)
                    yield from yield_locations_with_path(v, path)
                    path.pop()
        elif k == 'location':
            path.append('location')
            yield model[k], path
            path.pop()


def yield_points_from_location(location, path=[]):
    print(location)
    if location["type"] == "range":
        path.append('start')
        yield from yield_points_from_location(location["start"], path)
        path.pop()
        path.append('end')
        yield from yield_points_from_location(location["end"], path)
        path.pop()
    elif location["type"] == "point":
        yield location, path


def identify_unsorted_locations(model):
    """
    Identify first 2 locations that are not sorted from the model.
    Prerequisites:
    - no uncertainties;
    - greater or equal than 0 positions.

    :param model: Description model.
    :return:
    """
    current_point = 0
    current_point_path = []
    for location, location_path in yield_locations_with_path(model):
        if location.get('uncertain'):
            raise Exception("Uncertain point encountered.")
        for point, point_path in yield_points_from_location(location):
            if point.get('uncertain'):
                raise Exception("Uncertain point encountered.")
            else:
                if current_point > point['position']:
                    return current_point_path, location_path + point_path
                else:
                    current_point_path = location_path + point_path
    return None, None


def are_locations_sorted(model):
    """
    Checking if the locations present in the model are sorted.
    Prerequisites:
    - no uncertainties;
    - greater or equal than 0 positions.

    :param model: Description model.
    :return: True if the locations are sorted, False otherwise.
    """
    current_point = 0
    for location, location_path in yield_locations_with_path(model):
        for point, point_path in yield_points_from_location(location):
            if point.get('uncertain'):
                raise Exception("Uncertain point encountered.")
            else:
                if current_point > point["position"]:
                    return False
                else:
                    current_point = point["position"]
    return True


def yield_all_locations_models(model_1, model_2):
    """
    Yields all the location models present in the models, at all their levels.
    For (10_?)_30 there will be models for (10_?)_30, (10_?), 10, ?, and 30.
    The input models should have the same

    :param model_1: description model
    :param model_2: description model
    """
    for k in model_1.keys():
        if k in ['variants', 'deleted', 'inserted']:
            if isinstance(model_1[k], list):
                for i, v in enumerate(model_1[k]):
                    yield from yield_all_locations_models(model_1[k][i], model_2[k][i])
            elif isinstance(model_1[k], dict):
                yield from yield_all_locations_models(model_1[k], model_2[k])
        elif k == 'location':
            yield model_1[k], model_2[k]
            if model_1[k]["type"] == "range":
                yield model_1[k]["start"], model_2[k]["start"]
                if model_1[k]["start"]["type"] == "range":
                    yield model_1[k]["start"]["start"], model_2[k]["start"]["start"]
                    yield model_1[k]["start"]["end"], model_2[k]["start"]["end"]
                if model_1[k]["end"]["type"] == "range":
                    yield model_1[k]["end"]["start"], model_2[k]["end"]["start"]
                    yield model_1[k]["end"]["end"], model_2[k]["end"]["end"]
                yield model_1[k]["end"], model_2[k]["end"]


def contains_uncertain_locations(model):
    """
    Goes through model locations to see if any is uncertain.

    :param model: description model
    :return: True when the first uncertain location if encountered
             and False if none is encountered.
    """
    for location in yield_all_locations(model):
        if location.get('uncertain'):
            return True
    return False


def check_overlap(variants):
    sorted_variants = sort_variants(variants)


def check_points_2():
    pass


def check_locations(hgvs_model, internal_model, references=None):
    reference = get_reference(hgvs_model)
    for location_to_check, location_to_report in yield_all_locations_models(hgvs_model, internal_model):
        print('--')
        print(location_to_check)
        print(location_to_report)
    # check_points(hgvs_model, internal_model, reference)

