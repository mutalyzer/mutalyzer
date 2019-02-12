from schema import Schema, And, Or, Optional

location_offset_exact = Schema({'value': int})

location_offset_uncertain = Schema({'uncertain': And(bool, True)})

location_offset = Schema(Or(location_offset_exact, location_offset_uncertain))

location_point_exact = Schema({'type':                  And(str, 'point'),
                               'position':              int,
                               Optional('offset'):      location_offset,
                               Optional('outside_cds'): And(str, Or('downstream', 'upstream'))})

location_point_uncertain = Schema({'type':                  And(str, 'point'),
                                   'uncertain':             And(bool, True),
                                   Optional('offset'):      location_offset,
                                   Optional('outside_cds'): And(str, Or('downstream', 'upstream'))})

location_point = Schema(Or(location_point_exact, location_point_uncertain))

location_range_uncertain = Schema({'type':      And(str, 'range'),
                                   'uncertain': And(bool, True),
                                   'start':     location_point,
                                   'end':       location_point})

location_range = Schema({'type':  And(str, 'range'),
                         'start': Or(location_point, location_range_uncertain),
                         'end':   Or(location_point, location_range_uncertain)})

location = Schema(Or(location_point, location_range))

insertion_location = Schema({'source':             And(str, Or('observed', 'reference')),
                             'location':           location,
                             Optional('inverted'): bool})

insertion_sequence = Schema({'source':   And(str, 'description'),
                             'sequence': str})

insertion = Schema(Or(insertion_location, insertion_sequence))

deletion_sequence = Schema({'source': And(str, 'description'),
                            'sequence': str})

deletion_length = Schema({'source': And(str, 'description'),
                          'length': int})

deletion = Schema(Or(deletion_sequence, deletion_length))

variant = Schema({'type':     And(str, Or('equal', 'delins', 'inv', 'substitution', 'del', 'ins', 'dup', 'con')),
                  'location': location,
                  Optional('insertions'): [insertion],
                  Optional('deleted'):    [deletion]})

variants = Schema(Or([{'type': And(str, 'equal')}], [variant]))
