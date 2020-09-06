import copy

from normalizer.util import get_end, set_start


def update_interval(range_location, insertion=False):
    if (
        insertion
        and range_location["start"]["type"] == "point"
        and not range_location["start"].get("uncertain")
    ):
        range_location["start"]["position"] += 1
    elif (
        not range_location["end"].get("uncertain")
        and range_location["end"]["type"] == "point"
    ):
        range_location["end"]["position"] += 1


def point_to_range(point_location):
    """
    Convert a point location to a range location.

    :param point_location: A point location model complying object.
    :return: A range location model complying object.
    """
    if point_location.get("uncertain"):
        return {
            "type": "range",
            "start": {"type": "point", "uncertain": True},
            "end": {"type": "point", "uncertain": True},
        }
    return {
        "type": "range",
        "start": {"type": "point", "position": point_location["position"]},
        "end": {"type": "point", "position": point_location["position"] + 1},
    }


def location_to_internal_indexing(location, insertion=False):
    """
    """
    if location["type"] == "range":
        new_location = copy.deepcopy(location)
        update_interval(new_location, insertion)
    elif location["type"] == "point":
        new_location = point_to_range(location)
    else:
        # Should never happen. TODO: Maybe raise an error?
        pass
    return new_location


def variant_to_internal_indexing(variant):
    new_variant = copy.deepcopy(variant)
    if variant['type'] == "insertion":
        new_variant["location"] = location_to_internal_indexing(
            location=new_variant["location"], insertion=True)
    else:
        new_variant["location"] = location_to_internal_indexing(
            location=new_variant["location"])
    if new_variant.get("deleted"):
        for deleted in new_variant["deleted"]:
            if deleted.get("location"):
                deleted["location"] = location_to_internal_indexing(
                    deleted["location"])
    if new_variant.get("inserted"):
        for inserted in new_variant["inserted"]:
            if inserted.get("location"):
                inserted["location"] = location_to_internal_indexing(
                    inserted["location"])
    return new_variant


def to_internal_indexing(variants):
    """
    Convert the variants list to its deletion insertion only
    equivalent. It considers that internal indexing is employed.
    """
    new_variants = []
    for variant in variants:
        new_variants.append(variant_to_internal_indexing(variant))

    return new_variants
