import json

from mutalyzer.util import add_to_dict, create_exact_range_model


def generate_reference_selector_model(selector_model):
    reference_selector_model = {
        "id": selector_model["id"],
        "type": selector_model["type"],
        "location": create_exact_range_model(
            selector_model["exon"][0][0], selector_model["exon"][-1][-1]
        ),
    }
    add_to_dict(reference_selector_model, selector_model, "inverted")
    features = []
    for i, exon in enumerate(selector_model["exon"], start=1):
        features.append(
            {"type": "exon", "id": str(i), "location": create_exact_range_model(*exon)}
        )
    if selector_model.get("cds"):
        features.append(
            {
                "type": "cds",
                "id": str(i),
                "location": create_exact_range_model(*selector_model["cds"]),
            }
        )
    reference_selector_model["features"] = features

    return reference_selector_model


def generate_gene_model(selector_model):
    return [
        {
            "type": "gene",
            "id": "g1",
            "features": [generate_reference_selector_model(selector_model)],
        }
    ]


def generate_references(selector_model):
    record = {
        "annotations": {
            "type": "record",
            "id": "R1",
            "qualifiers": {"mol_type": "dna"},
            "features": generate_gene_model(selector_model),
        },
    }

    return {"R1": record, "reference": record}


def append_transcript(record, selector_model):
    for r in record:
        record[r]["annotations"]["features"][0]["features"].append(
            generate_reference_selector_model(selector_model)
        )


if __name__ == "__main__":
    pass
