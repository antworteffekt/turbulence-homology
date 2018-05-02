
import pytest
import pickle
import config


def load(s, var_name):

    ci = config.LES
    program_out = ci.PROGRAM_OUT
    test_out = ci.TEST_OUT

    fname_prog = '%s/%s/%s_structure_sizes_01.pickle' % (program_out, s, var_name)
    with open(fname_prog, 'r') as f:
        prog = pickle.load(f)
    fname_test = '%s/%s_LES_%s.pickle' % (test_out, var_name, s)
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

def test_SP1_z_wind():

    prog_data, test_data = load('SP1', 'z_wind')
    compare(prog_data, test_data)


def test_SP2_z_wind():

    prog_data, test_data = load('SP2', 'z_wind')
    compare(prog_data, test_data)

def test_SP3_z_wind():

    prog_data, test_data = load('SP3', 'z_wind')
    compare(prog_data, test_data)

def test_SP4_z_wind():

    prog_data, test_data = load('SP4', 'z_wind')
    compare(prog_data, test_data)

###########################################################################
# Second test block: y_wind

def test_SP1_y_wind():

    prog_data, test_data = load('SP1', 'y_wind')
    compare(prog_data, test_data)


def test_SP2_y_wind():

    prog_data, test_data = load('SP2', 'y_wind')
    compare(prog_data, test_data)

def test_SP3_y_wind():

    prog_data, test_data = load('SP3', 'y_wind')
    compare(prog_data, test_data)

def test_SP4_y_wind():

    prog_data, test_data = load('SP4', 'y_wind')
    compare(prog_data, test_data)

###########################################################################
# Third test block: x_wind

def test_SP1_x_wind():

    prog_data, test_data = load('SP1', 'x_wind')
    compare(prog_data, test_data)


def test_SP2_x_wind():

    prog_data, test_data = load('SP2', 'x_wind')
    compare(prog_data, test_data)

def test_SP3_x_wind():

    prog_data, test_data = load('SP3', 'x_wind')
    compare(prog_data, test_data)

def test_SP4_x_wind():

    prog_data, test_data = load('SP4', 'x_wind')
    compare(prog_data, test_data)