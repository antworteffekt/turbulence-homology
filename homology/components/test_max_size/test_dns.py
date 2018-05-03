
import pytest
import pickle
import config


def load(s, var_name):

    ci = config.DNS
    program_out = ci.PROGRAM_OUT
    test_out = ci.TEST_OUT

    fname_prog = '%s/%s/%s_max_structure_sizes.pickle' % (program_out, s, var_name)
    with open(fname_prog, 'r') as f:
        prog = pickle.load(f)
    fname_test = '%s/max_structure_DNS_%s.pickle' % (test_out, s)
    with open(fname_test, 'r') as f:
        test = pickle.load(f)

    return (prog, test)


def compare(prog_data, test_data):

    assert prog_data.keys() == test_data.keys()

    for d_1, d_2 in zip(prog_data.values(), test_data.values()):
        for k, v in d_1.items():
            assert d_1[k].keys() == d_2[k].keys()
            for d in d_1[k].keys():
                assert len(d_1[k][d]) == len(d_2[k][d])
                assert set(d_1[k][d]) == set(d_2[k][d])


def test_small():

    prog_data, test_data = load('SMALL', 'z_wind')
    compare(prog_data, test_data)


def test_medium():

    prog_data, test_data = load('MEDIUM', 'z_wind')
    compare(prog_data, test_data)