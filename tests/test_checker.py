import pytest
from normalizer.checker import is_overlap, is_sorted

VARIANTS = {'4_4delinsAA': {'type': 'deletion_insertion',
                            'source': 'reference',
                            'location': {'type': 'range',
                                         'start': {'type': 'point',
                                                   'position': 4},
                                         'end': {'type': 'point',
                                                 'position': 4}},
                            'inserted': [{'source': 'description',
                                          'sequence': 'AA'}]},
            '10del': {'type': 'deletion',
                      'source': 'reference',
                      'location': {
                          'type': 'point',
                          'position': 10}},
            '10_20=': {'type': 'equal',
                       'source': 'reference',
                       'location': {'type': 'range',
                                    'start': {'type': 'point',
                                              'position': 10},
                                    'end': {'type': 'point',
                                            'position': 20}}},
            '(10_15)del': {'type': 'deletion',
                           'source': 'reference',
                           'location': {'type': 'range',
                                        'uncertain': True,
                                        'start': {
                                            'type': 'point',
                                            'position': 10
                                        },
                                        'end': {'type': 'point',
                                                'position': 15}}},
            '(10_20)_(30_40)del': {'type': 'deletion',
                                   'source': 'reference',
                                   'location': {'type': 'range',
                                                'start': {'type': 'range',
                                                          'uncertain': True,
                                                          'start': {'type': 'point',
                                                                    'position': 10},
                                                          'end': {'type': 'point',
                                                                  'position': 20}},
                                                'end': {'type': 'range',
                                                        'uncertain': True,
                                                        'start': {'type': 'point',
                                                                  'position': 30},
                                                        'end': {'type': 'point',
                                                                'position': 40}}}},
            '(?_10)_(100_?)del': {'type': 'deletion',
                                  'source': 'reference',
                                  'location': {'type': 'range',
                                               'start': {'type': 'range',
                                                         'uncertain': True,
                                                         'start': {'type': 'point',
                                                                   'uncertain': True},
                                                         'end': {'type': 'point',
                                                                 'position': 10}},
                                               'end': {'type': 'range',
                                                       'uncertain': True,
                                                       'start': {'type': 'point',
                                                                 'position': 100},
                                                       'end': {'type': 'point',
                                                               'uncertain': True}}}}
            }


IS_SORTED = [(['4_4delinsAA', '10del'], True),
             (['10del', '(10_20)_(30_40)del'], True),
             (['10del', '4_4delinsAA'], False),
             (['10del', '4_4delinsAA'], False)]

IS_OVERLAP = [(['4_4delinsAA', '10del'], False),
              (['(10_15)del', '10del'], True),
              (['10_20=', '10del'], True),]


def get_variants(tests):
    output = []
    for test in tests:
        test_list = []
        for variant in test[0]:
            test_list.append(VARIANTS[variant])
        output.append((test_list, test[1]))
    return output


@pytest.mark.parametrize('variants, expected_result',
                         get_variants(IS_SORTED))
def test_is_sorted(variants, expected_result):
    assert is_sorted(variants) == expected_result


@pytest.mark.parametrize('variants, expected_result',
                         get_variants(IS_OVERLAP))
def test_is_overlap(variants, expected_result):
    assert is_overlap(variants) == expected_result


def test_semantic_check():
    pass


def test_reference_check():
    pass
