import configparser
import os
from pathlib import Path

from _collections import OrderedDict
from mutalyzer_mutator.util import reverse_complement


def create_exact_point_model(point):
    return {"type": "point", "position": point}


def create_exact_range_model(start, end):
    return {
        "type": "range",
        "start": create_exact_point_model(start),
        "end": create_exact_point_model(end),
    }


def get_start(model):
    """
    Get the start position of a (feature) location. For point locations
    the position value is returned. In case of uncertain start end,
    the minimum is returned.
    """
    if model.get("location"):
        model = model["location"]
    if model["type"] == "range":
        if model["start"].get("uncertain"):
            return get_start(model["start"])
        else:
            return model["start"]["position"]
    elif model["type"] == "point":
        return model["position"]


def set_start(location, position):
    location["start"]["position"] = position


def get_end(model):
    """
    Get the end position of a (feature) location. For point locations
    the position value is returned. In case of uncertain end range,
    the maximum is returned.
    """
    if model.get("location"):
        model = model["location"]
    if model["type"] == "range":
        return model["end"]["position"]
    elif model["type"] == "point":
        return model["position"]


def set_end(location, position):
    location["end"]["position"] = position


def get_location_length(location):
    return abs(get_end(location) - get_start(location))


def get_inserted_length(inserted):
    length = 0
    for insert in inserted:
        length += get_location_length(insert["location"])
    return length


def get_location_as_list(location):
    return [get_start(location), get_end(location)]


def sort_variants(variants):
    return sorted(
        variants,
        key=lambda variant: (
            get_start(variant["location"]),
            get_end(variant["location"]),
        ),
    )


def roll(s, first, last):
    """
    Determine the variability of a variant by looking at cyclic
    permutations. Not all cyclic permutations are tested at each time, it
    is sufficient to check ``aW'' if ``Wa'' matches (with ``a'' a letter,
    ``W'' a word) when rolling to the left for example.
        >>> roll('abbabbabbabb', 4, 6)
        (3, 6)
        >>> roll('abbabbabbabb', 5, 5)
        (0, 1)
        >>> roll('abcccccde', 4, 4)
        (1, 3)
    @arg s: A reference sequence.
    @type s: any sequence type
    @arg first: First position of the pattern in the reference sequence.
    @type first: int
    @arg last: Last position of the pattern in the reference sequence.
    @type last: int
    @return: tuple:
        - left  ; Amount of positions that the pattern can be shifted to
                  the left.
        - right ; Amount of positions that the pattern can be shifted to
                  the right.
    @rtype: tuple(int, int)
    """
    pattern = s[first - 1 : last]  # Extract the pattern
    pattern_length = len(pattern)

    # Keep rolling to the left as long as a cyclic permutation matches.
    minimum = first - 2
    j = pattern_length - 1
    while minimum > -1 and s[minimum] == pattern[j % pattern_length]:
        j -= 1
        minimum -= 1

    # Keep rolling to the right as long as a cyclic permutation matches.
    maximum = last
    j = 0
    while maximum < len(s) and s[maximum] == pattern[j % pattern_length]:
        j += 1
        maximum += 1

    return first - minimum - 2, maximum - last


def print_time_information(time_stamps):
    for i in range(1, len(time_stamps)):
        print(
            "{:<30}: {:2.6f}".format(
                time_stamps[i][0], time_stamps[i][1] - time_stamps[i - 1][1]
            )
        )
    print("{:<30}: {:2.6f}".format("TOTAL", time_stamps[-1][1] - time_stamps[0][1]))


def get_time_information(time_stamps):
    output = OrderedDict()
    for i in range(1, len(time_stamps)):
        output[time_stamps[i][0]] = "{:2.6f}".format(
            time_stamps[i][1] - time_stamps[i - 1][1]
        )
    output["TOTAL"] = "{:2.6f}".format(time_stamps[-1][1] - time_stamps[0][1])
    return output


def sort_location_tuples(locations):
    sorted_locations = sorted([i for j in locations for i in j])
    return list(zip(sorted_locations[0::2], sorted_locations[1::2]))


def string_k_v(width, key, value):
    return " {k:<{w}} : {v}\n".format(w=width, k=key, v=value)


def add_to_dict(d, source_d, k):
    if source_d.get(k):
        d[k] = source_d[k]


def add_msg(dictionary, message_type, message):
    if dictionary.get(message_type) is None:
        dictionary[message_type] = []
    dictionary[message_type].append(message)


def configuration():
    if (
        os.environ.get("MUTALYZER_SETTINGS")
        and Path(os.environ["MUTALYZER_SETTINGS"]).is_file()
    ):
        with open(os.environ["MUTALYZER_SETTINGS"]) as f:
            configuration_content = "[config]\n" + f.read()

        loaded_settings = configparser.ConfigParser()
        loaded_settings.optionxform = str
        loaded_settings.read_string(configuration_content)
        loaded_settings = {
            sect: dict(loaded_settings.items(sect))
            for sect in loaded_settings.sections()
        }["config"]

        return loaded_settings


def log_dir():
    settings = configuration()
    if settings and settings.get("MUTALYZER_LOG_DIR"):
        return eval(settings["MUTALYZER_LOG_DIR"])
    else:
        return "/tmp/mutalyzer.log"


def cache_dir():
    settings = configuration()
    if settings and settings.get("MUTALYZER_CACHE_DIR"):
        return eval(settings["MUTALYZER_CACHE_DIR"])


def set_by_path(dictionary, path, value):
    nested_dictionary = dictionary
    for k in path[:-1]:
        nested_dictionary = nested_dictionary[k]
    nested_dictionary[path[-1]] = value


def updated_by_path(dictionary, path, value):
    nested_dictionary = dictionary
    for k in path[:-1]:
        nested_dictionary = nested_dictionary[k]
    nested_dictionary[path[-1]].update(value)


def get_submodel_by_path(dictionary, path):
    nested_dictionary = dictionary
    for k in path:
        nested_dictionary = nested_dictionary[k]
    return nested_dictionary


def reverse_path(model, path):
    new_path = []
    for i, p in enumerate(path):
        if isinstance(p, int):
            new_path.append(len(get_submodel_by_path(model, path[:i])) - 1 - p)
        elif p == "start":
            new_path.append("end")
        elif p == "end":
            new_path.append("start")
        else:
            new_path.append(p)
    return new_path


def check_errors(fn):
    def wrapper(self):
        if not self.errors:
            return fn(self)
        if self.errors and self.stop_on_errors:
            raise Exception(str(self.errors))

    return wrapper


def slice_sequence(location, sequence):
    return sequence[get_start(location) : get_end(location)]


def construct_sequence(slices, sequences):
    seq = ""
    for slice in slices:
        if slice.get("sequence"):
            slice_seq = slice["sequence"]
        elif slice.get("location"):
            if isinstance(slice["source"], str):
                source = slice["source"]
            elif isinstance(slice["source"], dict):
                source = slice["source"]["id"]
            slice_seq = slice_sequence(slice["location"], sequences[source])
        else:
            raise Exception("Unrecognized slice", slice)
        if slice.get("repeat_number") and slice["repeat_number"].get("type") == "point":
            slice_seq = slice_seq * slice["repeat_number"]["value"]
        if slice.get("inverted"):
            slice_seq = reverse_complement(slice_seq)
        seq += slice_seq

    return seq


def get_inserted_sequence(variant, sequences):
    if variant.get("inserted"):
        return construct_sequence(variant["inserted"], sequences)
    return ""


def is_dna(sequence):
    for s in sequence:
        if s not in "ATCGBDHKMNRSVWY":
            return False

    return True


def is_rna(sequence):
    for s in sequence:
        if s not in "acgubdhkmnrsvwy":
            return False

    return True


def point_in_insertion(model, path):
    if (
        len(path) == 4
        and path[0] == "variants"
        and model["variants"][path[1]].get("type") == "insertion"
        and path[2] == "location"
    ):
        return True
    return False


def slice_seq(seq, slices):
    output = ""
    for s in slices:
        output += seq[s[0] : s[1]]
    return output

