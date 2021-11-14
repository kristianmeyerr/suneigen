#! python

import suneigen4py
import numpy as np

def test_cvodes_version():
    assert suneigen4py.__version__ == "0.0.1"

def ode(t, x):
    rhs = np.array([5.0, 6.0])
    return rhs

def ode_x0(t):
    x0 = np.array([3.0, 4.0])
    return x0

def test_cvodes_inputs():

    assert(suneigen4py.add(5,5)==12)
    print("test")


