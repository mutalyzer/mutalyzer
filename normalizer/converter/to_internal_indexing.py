import copy


def update_range_points(range_location, insertion=False):
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
    if range_location["start"]["type"] == "range":
        update_range_points(range_location["start"], insertion)
    if range_location["end"]["type"] == "range":
        update_range_points(range_location["end"], insertion)


def point_to_range(point_location):
    """
    Convert a point location to a range location.

    :param point_location: A point location model complying object.
    :return: A range location model complying object.
    """
    if point_location.get("uncertain"):
        location = {
            "type": "range",
            "start": {"type": "point", "uncertain": True},
            "end": {"type": "point", "uncertain": True},
        }
    else:
        location = {
            "type": "range",
            "start": {"type": "point", "position": point_location["position"]},
            "end": {"type": "point", "position": point_location["position"] + 1},
        }
    if point_location.get("shift"):
        location["start"]["shift"] = point_location["shift"]
        location["end"]["shift"] = point_location["shift"]
    return location


def location_to_internal_indexing(location, insertion=False):
    if location["type"] == "range":
        new_location = copy.deepcopy(location)
        update_range_points(new_location, insertion)
    elif location["type"] == "point":
        new_location = point_to_range(location)
    else:
        # Should never happen. TODO: Maybe raise an error?
        pass
    return new_location


def variant_to_internal_indexing(variant):
    new_variant = copy.deepcopy(variant)
    if variant.get("type") == "insertion":
        new_variant["location"] = location_to_internal_indexing(
            location=new_variant["location"], insertion=True
        )
    else:
        new_variant["location"] = location_to_internal_indexing(
            location=new_variant["location"]
        )
    if new_variant.get("deleted"):
        for deleted in new_variant["deleted"]:
            if deleted.get("location"):
                deleted["location"] = location_to_internal_indexing(deleted["location"])
    if new_variant.get("inserted"):
        for inserted in new_variant["inserted"]:
            if inserted.get("location"):
                inserted["location"] = location_to_internal_indexing(
                    inserted["location"]
                )
            if inserted.get("coordinate_system"):
                inserted["coordinate_system"] = "i"
    return new_variant


def variants_to_internal_indexing(variants):
    new_variants = []
    for variant in variants:
        new_variants.append(variant_to_internal_indexing(variant))

    return new_variants


def to_internal_indexing(model):
    new_model = copy.deepcopy(model)
    new_model["coordinate_system"] = "i"

    if new_model.get("variants"):
        new_model["variants"] = variants_to_internal_indexing(model["variants"])

    return new_model
