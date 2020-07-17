import json

from schema import And, Optional, Or, Schema, Use


def is_valid_IUPAC_NT(sequence):
    for nt in sequence:
        if nt not in (
            "A",
            "C",
            "G",
            "T",
            "U",
            "Y",
            "K",
            "M",
            "S",
            "W",
            "B",
            "D",
            "H",
            "V",
            "N",
        ):
            return False
    return True


location_offset = Schema(Or({"value": int}, {"uncertain": And(bool, True)}))

location_point_exact = Schema(
    {
        "type": And(str, "point"),
        "position": And(int, lambda n: n >= 0),
        Optional("offset"): location_offset,
        Optional("outside_cds"): And(str, Or("downstream", "upstream")),
    }
)

location_point_uncertain = Schema(
    {
        "type": And(str, "point"),
        "uncertain": And(bool, True),
        Optional("offset"): location_offset,
        Optional("outside_cds"): And(str, Or("downstream", "upstream")),
    }
)

location_point = Schema(Or(location_point_exact, location_point_uncertain))

location_range_uncertain = Schema(
    {
        "type": And(str, "range"),
        "uncertain": And(bool, True),
        "start": location_point,
        "end": location_point,
    }
)

location_range = Schema(
    {
        "type": And(str, "range"),
        "start": Or(location_point, location_range_uncertain),
        "end": Or(location_point, location_range_uncertain),
    }
)

location = Schema(Or(location_point, location_range))

insertion_location = Schema(
    {
        "source": And(str, len),
        "location": location,
        Optional("inverted"): And(bool, True),
    }
)

sequence = Schema(
    {
        "source": And(str, "description"),
        "sequence": And(str, len, Use(str.upper), lambda s: is_valid_IUPAC_NT(s)),
    }
)

insertion = Schema(Or(insertion_location, sequence))

deletion_length = Schema(
    {"source": And(str, "description"), "length": And(int, lambda n: n > 0)}
)

deletion = Schema(Or(sequence, deletion_length))

variant = Schema(
    {
        "type": And(
            str,
            Or(
                "equal",
                "deletion_insertion",
                "inversion",
                "substitution",
                "deletion",
                "insertion",
                "duplication",
                "conversion",
            ),
        ),
        "location": location,
        "source": And(str, "reference"),
        Optional("inserted"): [insertion],
        Optional("deleted"): [deletion],
    }
)

variants = Schema(
    Or([{"type": And(str, "equal"), "source": And(str, "reference")}], [variant])
)


print(json.dumps(variant.json_schema(""), indent=2))

s = Schema({"test": str, "nested": {Optional("other"): str}})
# print(json.dumps(s.json_schema("https://example.com/my-schema.json")))
# print(json.dumps(s.json_schema(), indent=2))
