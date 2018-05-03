
import pytest
import pickle
import config


def load(s, var_name):

    ci = config.DALES
    program_out = ci.PROGRAM_OUT
    test_out = ci.TEST_OUT

    fname_prog = '%s/%s_imicro2/%s_structure_sizes_01.pickle' % (program_out, s, var_name)
    with open(fname_prog, 'r') as f:
        prog = pickle.load(f)
    fname_test = '%s/%s_DALES_%s.pickle' % (test_out, var_name, s)
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


###########################################################################
# First test block: z_wind

def test_20130401_z_wind():

    prog_data, test_data = load('20130401', 'z_wind')
    compare(prog_data, test_data)

def test_20130716_z_wind():

    prog_data, test_data = load('20130716', 'z_wind')
    compare(prog_data, test_data)

def test_20130717_z_wind():

    prog_data, test_data = load('20130717', 'z_wind')
    compare(prog_data, test_data)

def test_20140717_z_wind():

    prog_data, test_data = load('20140717', 'z_wind')
    compare(prog_data, test_data)

###########################################################################
# Second test block: y_wind

def test_20130401_y_wind():

    prog_data, test_data = load('20130401', 'y_wind')
    compare(prog_data, test_data)

def test_20130716_y_wind():

    prog_data, test_data = load('20130716', 'y_wind')
    compare(prog_data, test_data)

def test_20130717_y_wind():

    prog_data, test_data = load('20130717', 'y_wind')
    compare(prog_data, test_data)

def test_20140717_y_wind():

    prog_data, test_data = load('20140717', 'y_wind')
    compare(prog_data, test_data)

###########################################################################
# Third test block: x_wind
def test_20130401_x_wind():

    prog_data, test_data = load('20130401', 'x_wind')
    compare(prog_data, test_data)

def test_20130716_x_wind():

    prog_data, test_data = load('20130716', 'x_wind')
    compare(prog_data, test_data)

def test_20130717_x_wind():

    prog_data, test_data = load('20130717', 'x_wind')
    compare(prog_data, test_data)

def test_20140717_x_wind():

    prog_data, test_data = load('20140717', 'x_wind')
    compare(prog_data, test_data)