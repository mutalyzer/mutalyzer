import mutalyzer.errors as errors
import mutalyzer.infos as infos

from .converter.to_rna import get_position_type
from .description_model import yield_sub_model
from .util import construct_sequence, get_end, get_start, sort_variants


def are_sorted(variants):
    """
    Check if the provided variants list is sorted.
    """
    sorted_variants = sort_variants(variants)
    if sorted_variants == variants:
        return True
    else:
        False


def is_overlap(variants):
    """
    Check whether there is overlap between the variants.

    TODO: Add support for fuzzy (uncertain) locations.
    """
    sorted_variants = sort_variants(variants)
    positions = []
    for point, path in yield_sub_model(
        sorted_variants, ["location", "start", "end"], ["point"]
    ):
        if point.get("position") is not None:
            positions.append(point["position"])
    if positions:
        min_position = min(positions) - 1
        for variant in sorted_variants:
            if get_start(variant["location"]) <= min_position - 1:
                return True
            else:
                min_position = get_end(variant["location"])
    return False


def contains_uncertain_locations(model):
    """
    Goes through model locations to see if any is uncertain.

    :param model: description model
    :return: True when the first uncertain location if encountered
             and False if none is encountered.
    """
    for location, path in yield_sub_model(
        model, ["location", "start", "end"], ["point", "range"]
    ):
        if location.get("uncertain") or (
            location.get("offset") and location["offset"].get("uncertain")
        ):
            return True
    return False


def contains_insert_length(model):
    if model.get("variants"):
        for v in model["variants"]:
            if v.get("inserted"):
                if v["type"] not in ["duplication", "inversion"] or (
                    v["type"] in ["duplication", "inversion"] and len(v["inserted"]) > 1
                ):
                    for i in v["inserted"]:
                        if i.get("length"):
                            return True
    return False


def _in_adjacent_exons(start_i, end_i, exons):
    flatten_exons = [
        e for exon in exons[start_i[0] // 2 : end_i[0] // 2 + 1] for e in exon
    ]
    if len(flatten_exons) - len(set(flatten_exons)) == end_i[0] // 2 - start_i[0] // 2:
        return True
    else:
        return False


def splice_sites(variants, sequences, selector_model):
    def _starts_ends(v):
        start = v["location"]["start"]
        end = v["location"]["end"]
        if (
            start.get("shift")
            and start["shift"] > 0
            and end.get("shift")
            and start["shift"] == end["shift"]
        ):
            shift = start["shift"]
        else:
            shift = 0
        for idx in range(shift + 1):
            _start = start["position"] - idx
            if start["position"] == end["position"]:
                _end = _start
            else:
                _end = end["position"] - 1 - idx
            yield get_position_type(_start, selector_model["exon"]), get_position_type(
                _end, selector_model["exon"]
            )

    errors_splice = []
    infos_splice = []
    for i, v in enumerate(variants):
        if v.get("location"):
            path = ["variants", i]
            for start_i, end_i in _starts_ends(v):
                if start_i[0] % 2 == 1 and end_i[0] % 2 == 1:
                    # both start and end are in exons
                    if end_i[0] != start_i[0] and not _in_adjacent_exons(
                        start_i, end_i, selector_model["exon"]
                    ):
                        if v.get("inserted") and construct_sequence(
                            v["inserted"], sequences
                        ):
                            # insertion - error
                            errors_splice.append(errors.splice_site(["variants", i]))
                            break
                        else:
                            infos_splice.append(infos.splice_site_removed(path))
                elif start_i[0] % 2 == 0 and end_i[0] % 2 == 0:
                    # both start and end are in introns
                    if start_i[0] == end_i[0] and start_i[1] == end_i[1] == 0:
                        # same intron discarded - warning
                        infos_splice.append(infos.variant_discarded(path))
                    elif start_i[0] != end_i[0] and start_i[1] == 0 and end_i[1] == 0:
                        if v.get("inserted") and construct_sequence(
                            v["inserted"], sequences
                        ):
                            errors_splice.append(errors.splice_site(path))
                            break
                        else:
                            infos_splice.append(infos.splice_site_removed(path))
                    else:
                        errors_splice.append(errors.splice_site(path))
                        break
                else:
                    errors_splice.append(errors.splice_site(path))
                    break
    return errors_splice, infos_splice
