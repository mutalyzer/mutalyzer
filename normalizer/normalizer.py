from .description import Description
import json


def normalize(description_to_normalize):
    description = Description(description_to_normalize)

    description.normalize()
    print(json.dumps(description.infos, indent=2))
    print(json.dumps(description.errors, indent=2))

    return description.normalized_description
