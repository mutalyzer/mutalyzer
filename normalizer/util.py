

def get_start(location):
    """
    Get the start position of a location. For point locations the position
    value is returned. In case of uncertain start end, the minimum is returned.
    """
    if location['type'] == 'range':
        if location['start'].get('uncertain'):
            return get_start(location['start'])
        else:
            return location['start']['position']
    elif location['type'] == 'point':
        return location['position']


def set_start(location, start):
    location['start']['position'] = start


def get_end(location):
    """
    Get the end position of a location. For point locations the position value
    is returned. In case of uncertain end range, the maximum is returned.
    """
    if location['type'] == 'range':
        return location['end']['position']
    elif location['type'] == 'point':
        return location['position']


def update_position(location, start_end,  value):
    if location.get('type') == range:
        location[start_end] = {'type': 'point',
                               'position': value}


def sort_variants(variants):
    return sorted(variants,
                  key=lambda variant: get_start(variant['location']))


def to_description(variant):
    pass