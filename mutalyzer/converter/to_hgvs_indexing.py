import copy


def update_range_points(range_location, insertion=False):
    if (
        insertion
        and range_location["start"]["type"] == "point"
        and not range_location["start"].get("uncertain")
    ):
        range_location["start"]["position"] -= 1
    elif (
        not range_location["end"].get("uncertain")
        and range_location["end"]["type"] == "point"
    ):
        range_location["end"]["position"] -= 1
    if range_location["start"]["type"] == "range":
        update_range_points(range_location["start"], insertion)
    if range_location["end"]["type"] == "range":
        update_range_points(range_location["end"], insertion)


def location_to_hgvs_indexing(location, insertion=False):
    """"""
    new_location = copy.deepcopy(location)
    if (
        not insertion
        and new_location["start"].get("position") is not None
        and new_location["end"].get("position") is not None
        and new_location["start"]["position"] + 1 == new_location["end"]["position"]
    ):
        new_location = new_location["start"]
    else:
        update_range_points(new_location, insertion)

    return new_location


def variant_to_internal_indexing(variant):
    new_variant = copy.deepcopy(variant)
    if variant.get("type") == "insertion":
        new_variant["location"] = location_to_hgvs_indexing(
            location=new_variant["location"], insertion=True
        )
    else:
        if new_variant.get("location"):
            new_variant["location"] = location_to_hgvs_indexing(
                location=new_variant["location"]
            )
    if new_variant.get("deleted"):
        for deleted in new_variant["deleted"]:
            if deleted.get("location"):
                deleted["location"] = location_to_hgvs_indexing(deleted["location"])
    if new_variant.get("inserted"):
        for inserted in new_variant["inserted"]:
            if inserted.get("location"):
                inserted["location"] = location_to_hgvs_indexing(inserted["location"])
    return new_variant


def variants_to_internal_indexing(variants):
    new_variants = []
    for variant in variants:
        new_variants.append(variant_to_internal_indexing(variant))

    return new_variants


def to_hgvs_indexing(model):
    new_model = copy.deepcopy(model)
    new_model["coordinate_system"] = "x"

    if new_model.get("variants"):
        new_model["variants"] = variants_to_internal_indexing(model["variants"])

    return new_model
