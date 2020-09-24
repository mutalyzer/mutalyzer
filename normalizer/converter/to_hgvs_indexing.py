import copy


def location_to_hgvs_indexing(location, insertion=False):
    """
    """
    if insertion:
        new_location = copy.deepcopy(location)
        new_location["start"]["position"] -= 1
    elif location["start"]["position"] + 1 == location["end"]["position"]:
        new_location = copy.deepcopy(location["start"])
    else:
        new_location = copy.deepcopy(location)
        new_location["end"]["position"] -= 1

    return new_location


def variant_to_internal_indexing(variant):
    new_variant = copy.deepcopy(variant)
    if variant.get("type") == "insertion":
        new_variant["location"] = location_to_hgvs_indexing(
            location=new_variant["location"], insertion=True
        )
    else:
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
