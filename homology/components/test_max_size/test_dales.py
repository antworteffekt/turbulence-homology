
import pytest
import pickle
import config


def load(s, var_name):

    ci = config.DALES
    program_out = ci.PROGRAM_OUT
    test_out = ci.TEST_OUT

    fname_prog = '%s/%s/%s_max_structure_sizes.pickle' % (program_out, s, var_name)
    with open(fname_prog, 'r') as f:
        prog = pickle.load(f)
    fname_test = '%s/max_structure_DALES_%s.pickle' % (test_out, s)
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


###########################################################################
# First test block: z_wind

def test_20130401_z_wind():

    prog_data, test_data = load('20130401_imicro2', 'z_wind')
    compare(prog_data, test_data)

def test_20130716_z_wind():

    prog_data, test_data = load('20130716_imicro2', 'z_wind')
    compare(prog_data, test_data)

def test_20130717_z_wind():

    prog_data, test_data = load('20130717_imicro2', 'z_wind')
    compare(prog_data, test_data)

def test_20140717_z_wind():

    prog_data, test_data = load('20140717_imicro2', 'z_wind')
    compare(prog_data, test_data)
