import json

from .util import get_location_as_list
from .description import variant_to_description


def to_be_visualized(de_variants, hgvs_variants):
    output = []

    j = 0
    for de_variant in de_variants:
        if de_variant["type"] == "equal":
            output.append(
                {
                    "location": get_location_as_list(de_variant["location"]),
                    "type": "equal",
                }
            )
        else:
            o = {
                "location": get_location_as_list(de_variant["location"]),
                "type": hgvs_variants[j]["type"],
                "repr": variant_to_description(hgvs_variants[j]),
                "hover": False,
            }
            if hgvs_variants[j].get("inserted"):
                s = ""
                for inserted in hgvs_variants[j]["inserted"]:
                    if inserted.get("sequence"):
                        s += inserted.get("sequence")
                o["inserted"] = s
            output.append(o)
            j += 1
    return output
