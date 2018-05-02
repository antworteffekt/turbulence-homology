
import pytest
import pickle
import config


def load(s, var_name):

    ci = config.DNS
    program_out = ci.PROGRAM_OUT
    test_out = ci.TEST_OUT

    fname_prog = '%s/%s/%s_structure_sizes_01.pickle' % (program_out, s, var_name)
    with open(fname_prog, 'r') as f:
        prog = pickle.load(f)
    fname_test = '%s/%s_DNS_%s.pickle' % (test_out, var_name, s)
    with open(fname_test, 'r') as f:
        test = pickle.load(f)

    return (prog, test)


def compare(prog_data, test_data):

    assert prog_data.keys() == test_data.keys()

    for k, v in prog_data.items():
        assert prog_data[k].keys() == test_data[k].keys()
        for d in prog_data[k].keys():
            assert len(prog_data[k][d]) == len(test_data[k][d])
            assert set(prog_data[k][d]) == set(test_data[k][d])


def test_small():

    prog_data, test_data = load('SMALL', 'z_wind')
    compare(prog_data, test_data)


def test_medium():

    prog_data, test_data = load('MEDIUM', 'z_wind')
    compare(prog_data, test_data)