import bisect
from copy import deepcopy

from Bio.Seq import Seq
from mutalyzer_crossmapper import Coding, Genomic, NonCoding
from mutalyzer_mutator.util import reverse_complement

from ..description_model import (
    variant_to_description,
    variants_to_description,
    yield_sub_model,
)
from ..reference import (
    extract_feature_model,
    get_internal_selector_model,
    slice_to_selector,
    yield_locations,
)
from ..util import (
    construct_sequence,
    get_end,
    get_inserted_sequence,
    get_start,
    set_by_path,
    set_end,
    set_start,
)
from .to_hgvs_coordinates import genomic_to_point, reverse_strand_shift


def to_rna_reference_model(reference_model, selector_id, transcribe=True):
    """
    Get the RNA reference model of the provided selector.

    1. Extract the tree corresponding to the selector from the model (including
    the parents).
    2. Slice the sequence.
    3. Update the model features locations using the crossmapper.

    TODO: Make sure everything is on the plus strand?

    :arg dict reference_model: Reference model.
    :arg str selector_id: Selector ID.
    :arg bool transcribe: Transcribe the sequence to RNA.
    :returns: RNA reference model.
    :rtype: dict
    """
    rna_model = {
        "annotations": deepcopy(
            extract_feature_model(reference_model["annotations"], selector_id)[0]
        ),
        "sequence": {
            "seq": str(
                Seq(slice_to_selector(reference_model, selector_id)).transcribe()
            ).lower()
            if transcribe
            else slice_to_selector(reference_model, selector_id)
        },
    }
    s_m = get_internal_selector_model(rna_model["annotations"], selector_id, True)
    x = NonCoding(s_m["exon"]).coordinate_to_noncoding

    new_start = x(s_m["exon"][0][0])[0] - 1
    new_end = x(s_m["exon"][-1][-1])[0]
    for location, f_type in yield_locations(rna_model["annotations"]):
        if f_type == "CDS":
            set_start(location, x(get_start(location))[0] - 1)
            set_end(location, x(get_end(location))[0] - 1)
        elif f_type == "exon":
            set_start(location, x(get_start(location))[0] - 1)
            set_end(location, x(get_end(location))[0] + x(get_end(location))[1] - 1)
        else:
            set_start(location, new_start)
            set_end(location, new_end)
    return rna_model


def get_position_type(position, exons, len_ss=2, len_as=5):
    """
    Get the position location within the exons/introns. Even numbers for
    introns and odd numbers for exons are returned. Empty introns are
    considered as well in the returned index. The second returned value
    represents a splice site (1, -1) or around a splice site (-2, 2) location,
    otherwise 0 (within an intron outside the splice (around) sites or
    within an exon).

    :arg int position: Zero-based position.
    :arg list exons: Zero-based half open exon positions list of tuples.
    :arg int len_ss: Splice site length.
    :arg int len_as: Around splice site length.
    :returns: Position type.
    :rtype: tuple
    """
    x = NonCoding(exons).coordinate_to_noncoding
    exons = _get_flatten_exons(exons)
    position_x = x(position)

    if position_x[1] == 0:
        return bisect.bisect_right(exons, position), 0
    elif 0 < abs(position_x[1]) <= len_ss:
        if position_x[1] > 0:
            return bisect.bisect_right(exons, position), 1
        else:
            return bisect.bisect_left(exons, position), -1
    elif len_ss < abs(position_x[1]) <= len_ss + len_as:
        if position_x[1] > 0:
            return bisect.bisect_right(exons, position), 2
        else:
            return bisect.bisect_left(exons, position), -2
    else:
        return bisect.bisect_left(exons, position), 0


def get_location_type(location, exons, len_ss=2, len_as=5):
    """
    Returns the location spanning with respect to the exons/introns. Currently
    the supported types are: same exon (start and end in the same exon),
    exon - exon (start and end in different exons), same intron,
    and intron - intron.

    :arg dict location: Location model.
    :arg list exons: Flatten exon positions.
    :returns: Location type within the exons/introns.
    :rtype: str
    """
    start_i = get_position_type(get_start(location), exons, len_ss, len_as)
    end_i = get_position_type(get_end(location) - 1, exons, len_ss, len_as)
    if get_start(location) == get_end(location):
        # this is an insertion
        if start_i[0] % 2 == 1:
            return "same exon"
        else:
            if start_i[1] == 0:
                return "same intron"
    elif start_i[0] % 2 == 1 and end_i[0] % 2 == 1:
        if start_i[0] == end_i[0]:
            return "same exon"
        else:
            return "exon exon"
    elif start_i[0] % 2 == 0 and end_i[0] % 2 == 0:
        if start_i[0] == end_i[0] and start_i[1] == 0:
            return "same intron"
        if start_i[0] != end_i[0] and start_i[1] == 0 and end_i[1] == 0:
            return "intron intron"
    elif start_i[0] % 2 == 1 and end_i[0] % 2 == 0:
        return "exon intron"
    elif start_i[0] % 2 == 0 and end_i[0] % 2 == 1:
        return "intron exon"


def _get_flatten_exons(exons):
    """
    Transform the exon list of tuples into a list of integers.

    :params list exons: Exons as a list of tuples.
    :return: Flattened exons list.
    :rtype: list
    """
    return [e for exon in exons for e in exon]


def _get_exon_start_position(position, exons):
    """
    Given an intronic position (start), get its appropriate exon position.

    :arg int position: Zero-based position.
    :arg list exons: Flattened exons list.
    :returns: Exon position.
    :rtype: int
    """
    return exons[bisect.bisect_right(exons, position)]


def _get_exon_end_position(position, exons):
    """
    Given an intronic position (end), get its appropriate exon position.

    :arg int position: Zero-based position.
    :arg list exons: Flattened exons list.
    :returns: Exon position.
    :rtype: int
    """
    return exons[bisect.bisect_left(exons, position) - 1]


def _set_start_to_exon(location, exons):
    """
    Update the location start position with its appropriate exon position.

    :arg dict location: Zero-based location model.
    :arg list exons: Flattened exons list.
    """
    set_start(location, _get_exon_start_position(get_start(location), exons))


def _set_end_to_exon(location, exons):
    """
    Update the location end position with its appropriate exon position.

    :arg dict location: Zero-based location model.
    :arg list exons: Flattened exons list.
    """
    set_end(location, _get_exon_end_position(get_end(location), exons))


def _trim_to_exons(variants, exons, sequences):
    """
    Update variants locations to the corresponding exons.
    Notes:
      - same intron locations are discarded;
      - splice sites checked should have been performed already.
    """
    new_variants = []
    for v in variants:
        new_v = deepcopy(v)
        if v.get("location"):
            location_type = get_location_type(v["location"], exons)
            if location_type == "intron intron" and not (
                v.get("inserted") and construct_sequence(v["inserted"], sequences)
            ):
                _set_start_to_exon(new_v["location"], _get_flatten_exons(exons))
                _set_end_to_exon(new_v["location"], _get_flatten_exons(exons))
                new_variants.append(new_v)
            elif location_type == "exon exon":
                new_variants.append(new_v)
            elif location_type == "same exon":
                new_variants.append(new_v)
    return new_variants


def to_rna_variants(variants, sequences, selector_model):
    """
    Convert coordinate delins variants to RNA.

    :arg list variants: Variants with coordinate locations.
    :arg list sequences: List with sequences dictionary.
    :arg dict selector_model: Selector model.
    :returns: Converted RNA variants.
    :rtype: dict
    """
    trimmed_variants = _trim_to_exons(variants, selector_model["exon"], sequences)

    x = NonCoding(selector_model["exon"]).coordinate_to_noncoding

    for variant in trimmed_variants:
        if variant.get("location"):
            set_start(variant["location"], x(get_start(variant))[0] - 1)
            set_end(
                variant["location"], x(get_end(variant))[0] + x(get_end(variant))[1] - 1
            )
            if variant.get("inserted"):
                variant["inserted"] = [
                    {
                        "source": "description",
                        "sequence": get_inserted_sequence(variant, sequences),
                    }
                ]
    return to_rna_sequences(trimmed_variants)


def to_rna_sequences(model):
    """
    Convert all the sequences present in the model to RNA.

    :args dict model: Description model.
    """
    for seq, path in yield_sub_model(model, ["sequence"]):
        set_by_path(model, path, str(Seq(seq).transcribe().lower()))
    return model


def _point_to_cds_coordinate(point, selector_model, crossmap):
    genomic_to_coordinate = Genomic().genomic_to_coordinate
    if selector_model.get("inverted"):
        if point.get("shift"):
            point["position"] -= point["shift"]
    coding = crossmap.coordinate_to_coding(point["position"], degenerate=True)
    if coding[2] == -1:
        return genomic_to_point(0)
    else:
        return genomic_to_point(genomic_to_coordinate(coding[0] + coding[1]))


def _get_inserted_sequence(insertion, sequences):
    if isinstance(insertion["source"], str):
        source = insertion["source"]
    elif isinstance(insertion["source"], dict):
        source = insertion["source"]["id"]
    return sequences[source][
        get_start(insertion["location"]) : get_end(insertion["location"])
    ]


def merge_inserted_to_string(inserted, sequences):
    inserted_value = ""
    for insertion in inserted:
        if insertion.get("sequence"):
            inserted_value += insertion.get("sequence")
        else:
            inserted_value += _get_inserted_sequence(insertion, sequences)
        if insertion.get("inverted"):
            inserted_value = reverse_complement(inserted_value)

    return {"source": "description", "sequence": inserted_value}


def variant_to_cds_coordinate(variant, sequences, selector_model, crossmap):
    new_variant = deepcopy(variant)

    location = new_variant["location"]

    if location["type"] == "range":
        location["start"] = _point_to_cds_coordinate(
            location["start"], selector_model, crossmap
        )
        location["end"] = _point_to_cds_coordinate(
            location["end"], selector_model, crossmap
        )
    else:
        location = _point_to_cds_coordinate(location, selector_model, crossmap)
    if new_variant.get("inserted"):
        new_variant["inserted"] = [
            {
                "source": "description",
                "sequence": get_inserted_sequence(variant, sequences),
            }
        ]
    new_variant["location"] = location
    return new_variant


def reverse_start_end(variants):
    for variant in variants:
        if variant.get("location") and variant["location"]["type"] == "range":
            location = variant["location"]
            location["start"], location["end"] = location["end"], location["start"]
            location["start"]["position"] -= 1
            location["end"]["position"] -= 1


def _get_cds_into_exons(exons, cds):
    l_index = bisect.bisect_right(exons, cds[0])
    r_index = bisect.bisect_left(exons, cds[1])
    return [cds[0]] + exons[l_index:r_index] + [cds[1]]


def _location_in_same_intron(location, exons):
    start_i = bisect.bisect_right(exons, get_start(location))
    end_i = bisect.bisect_left(exons, get_end(location))
    if start_i == end_i and start_i % 2 == 0:
        return True
    else:
        return False


def _splice_site_removal(location, exons):
    start_i = bisect.bisect_right(exons, get_start(location))
    end_i = bisect.bisect_left(exons, get_end(location))
    if end_i - start_i == 1:
        return True


def _get_exons_and_cds(selector_model):
    exons = [e for t in selector_model["exon"] for e in t]
    cds = [selector_model["cds"][0][0], selector_model["cds"][0][1]]
    if selector_model.get("inverted"):
        cds[0] = exons[0]
    else:
        cds[1] = exons[-1]
    return exons, cds


def _get_exons_and_cds_2(s_m):
    exons = [e for t in s_m["exon"] for e in t]
    cds = [s_m["cds"][0][0], s_m["cds"][0][1]]
    return exons, cds


def to_exon_positions(variants, exons, cds):
    exons = _get_cds_into_exons(exons, cds)
    new_variants = []
    for variant in variants:
        if (
            variant.get("type") == "deletion_insertion"
            and variant.get("location")
            and not _location_in_same_intron(variant["location"], exons)
            and not (get_start(variant) <= exons[0] and get_end(variant) <= exons[0])
        ):
            n_v = deepcopy(variant)
            exon_s = bisect.bisect(exons, get_start(n_v))
            if exon_s % 2 == 0 and exon_s < len(exons):
                n_v["location"]["start"]["position"] = exons[exon_s]

            exon_e = bisect.bisect(exons, get_end(n_v))
            if exon_e % 2 == 0 and exon_e < len(exons):
                n_v["location"]["end"]["position"] = exons[exon_e]

            new_variants.append(n_v)

    return new_variants


def _get_splice_site_hits(variants, exons, cds):
    hits = []
    for i, variant in enumerate(variants):
        if (
            variant.get("type") == "deletion_insertion"
            and variant.get("location")
            and _splice_site_removal(
                variant["location"], _get_cds_into_exons(exons, cds)
            )
        ):
            hits.append(i)
    return hits


def reverse_variants(variants, sequences):
    reversed_variants = deepcopy(variants)
    reverse_strand_shift(reversed_variants, sequences["reference"])
    reverse_start_end(reversed_variants)
    return reversed_variants


def to_rna_protein_coordinates(variants, sequences, selector_model):
    """
    Converts the locations to cds equivalent.

    :param variants: Variants with locations in the coordinate system.
    :param sequences: Sequences with their ids as keys.
    :param selector_model: Selector model according to which
                           the conversion is performed.
    """
    exons, cds = _get_exons_and_cds(selector_model)
    crossmap = Coding(selector_model["exon"], cds, selector_model["inverted"])

    if selector_model.get("inverted"):
        variants = reverse_variants(variants, sequences)

    splice_site_hits = _get_splice_site_hits(variants, exons, cds)

    coordinate_variants = to_exon_positions(variants, exons, cds)

    cds_variants = []
    for variant in coordinate_variants:
        cds_variants.append(
            variant_to_cds_coordinate(variant, sequences, selector_model, crossmap)
        )
    return cds_variants, splice_site_hits
