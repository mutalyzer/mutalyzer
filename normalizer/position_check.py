from .converter.to_internal_coordinates import get_reference
from .description import point_to_description


def check_location_sequence_boundary(hgvs_location, internal_location, reference):
    if hgvs_location['type'] == 'range':
        check_location_sequence_boundary(
            hgvs_location['start'], internal_location['start'], reference)
        check_location_sequence_boundary(
            hgvs_location['end'], internal_location['end'], reference)
    if hgvs_location ['type'] == 'point':
        if (internal_location['position'] > reference.get_length()
                or internal_location['position'] < 0):
            hgvs_location['errors'] = [{
                "code": "EOUTOFBOUNDARY",
                "details": "Position {} is out of sequence boundaries.".format(
                    point_to_description(hgvs_location))
            }]
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


def check_positions(hgvs_model, internal_model):
    reference = get_reference(hgvs_model)
    check_points(hgvs_model, internal_model, reference)

