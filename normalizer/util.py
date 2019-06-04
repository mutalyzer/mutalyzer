

def get_start(location):
    """
    Get the start position of a location. For point locations the position
    value is returned. In case of uncertain start end, the minimum is returned.
    """
    if location['type'] == 'range':
        if location['start'].get('uncertain'):
            return get_start(location['start'])
        else:
            return location['start']['position']
    elif location['type'] == 'point':
        return location['position']


def set_start(location, start):
    location['start']['position'] = start


def get_end(location):
    """
    Get the end position of a location. For point locations the position value
    is returned. In case of uncertain end range, the maximum is returned.
    """
    if location['type'] == 'range':
        return location['end']['position']
    elif location['type'] == 'point':
        return location['position']


def update_position(location, start_end,  value):
    if location.get('type') == range:
        location[start_end] = {'type': 'point',
                               'position': value}


def get_location_length(location):
    return abs(get_end(location) - get_start(location))


def get_inserted_length(inserted):
    length = 0
    for insert in inserted:
        length += get_location_length(insert['location'])
    return length


def sort_variants(variants):
    return sorted(variants,
                  key=lambda variant: get_start(variant['location']))


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
    pattern = s[first - 1:last]   # Extract the pattern
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
