from _collections import OrderedDict


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


def set_start(location, start):
    location["start"]["position"] = start


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


def update_position(location, start_end, value):
    if location.get("type") == range:
        location[start_end] = {"type": "point", "position": value}


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
    return sorted(variants, key=lambda variant: get_start(variant["location"]))


def to_description(variant):
    pass


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
