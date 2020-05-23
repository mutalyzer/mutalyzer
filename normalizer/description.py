

def get_selector_id(description_model):
    """
    Get the selector ID from the description model. At the moment, no nesting
    is supported.
    :param description_model: Provided by the HGVS description parser.
    :return: The ID of the selector, if provided, otherwise None.
    """
    if description_model.get('reference') and \
            description_model['reference'].get('selector') and \
            description_model['reference']['selector'].get('id'):
        return description_model['reference']['selector']['id']


def get_coordinate_system(description_model):
    if description_model.get('coordinate_system'):
        return description_model['coordinate_system']
