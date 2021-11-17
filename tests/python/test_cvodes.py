import suneigen4py
import numpy as np


def test_cvodes_version():
    assert suneigen4py.__version__ == "0.0.1"


def ode(t, x):
    rhs = np.array([0., 0., 0.])
    rhs[0] = -0.04 * x[0] + 1.0e4*x[1]*x[2]
    rhs[2] = 3.0e7 * x[1] * x[1]
    rhs[1] = -rhs[0] - rhs[2]
    return rhs


def jac(t, x):
    jac_ = np.array(
        [
            [-0.04,              1.0e4*x[2],    1.0e4*x[1]],
            [0.04, -1.0e4*x[2] - 6.0e7*x[1], -1.0e4 * x[1]],
            [0.0,              6.0e7 * x[1],           0.0]
        ]
    )
    return jac_


def test_cvodes_inputs():

    times = np.array([40.0, 400.0])
    x0 = np.array([1.0, 0.0, 0.0])
    options = {"AbsTol": 1e-12, "RelTol": 1e-08}

    t, x = suneigen4py.cvodes(ode, jac, x0, times, options)

    np.testing.assert_almost_equal(x[0,0], 7.158017e-01, decimal=4)
    np.testing.assert_almost_equal(x[0,1], 9.1850e-06, decimal=9)
    np.testing.assert_almost_equal(x[0,2], 2.8418e-01, decimal=4)
    np.testing.assert_almost_equal(x[1,0], 4.5053e-01, decimal=4)
    np.testing.assert_almost_equal(x[1,1], 3.2232e-06, decimal=9)
    np.testing.assert_almost_equal(x[1,2], 5.4946e-01, decimal=4)
