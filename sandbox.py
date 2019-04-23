from schema import Schema, And, Or, Optional, Use

feature_point = Schema({'position':            And(int, lambda n : n >= 0),
                        Optional('uncertain'): And(bool, True)})

feature = Schema({'features': Schema({}, ignore_extra_keys=True),
                  'start':    feature_point,
                  'end':      feature_point,
                  'type':     And(str, len)})


def validate(current):
    feature.validate(current)
    for next in current['features'].values():
        validate(next)


TEST = {'features': {'nested': {'features': {}, 'start': {'position': 0}, 'end': {'position': 0}, 'type': 'bla'}}, 'start': {'position': 0}, 'end': {'position': 0}, 'type': 'test'}

validate(TEST)
