import suneigen4py
import numpy as np
from scipy.sparse import csc_matrix

def test_cvodes_version():
    assert suneigen4py.__version__ == "0.0.1"


def ode(t, x, p, k, h):
    rhs = np.array([0., 0., 0.])
    rhs[0] = -0.04 * x[0] + 1.0e4*x[1]*x[2]
    rhs[2] = 3.0e7 * x[1] * x[1]
    rhs[1] = -rhs[0] - rhs[2]
    return rhs


def jac(t, x, p, k, h):
    jac_ = np.array(
        [
            [-0.04, 1.0e4*x[2], 1.0e4*x[1]],
            [0.04, (-1.0e4*x[2]) - (6.0e7*x[1]), -1.0e4*x[1]],
            [0.0, 6.0e7*x[1], 0.0]
        ]
    )
    jac_sp = csc_matrix(jac_)
    return jac_sp.indices, jac_sp.indptr, jac_sp.data


def test_cvodes_dense():

    times = np.array([40.0, 400.0])
    x0 = np.array([1.0, 0.0, 0.0])
    nnz = 7
    options = {"AbsTol": 1e-12, "RelTol": 1e-08, "linsol": "densE"}
    p = np.array([0.04, 1.0e4, 3.0e7])
    k = np.array([1.0])

    t, x = suneigen4py.cvodes(ode, jac, x0, times, nnz, p, k, options)

    np.testing.assert_almost_equal(x[0,0], 7.158017e-01, decimal=4)
    np.testing.assert_almost_equal(x[0,1], 9.1850e-06, decimal=9)
    np.testing.assert_almost_equal(x[0,2], 2.8418e-01, decimal=4)
    np.testing.assert_almost_equal(x[1,0], 4.5053e-01, decimal=4)
    np.testing.assert_almost_equal(x[1,1], 3.2232e-06, decimal=9)
    np.testing.assert_almost_equal(x[1,2], 5.4946e-01, decimal=4)


def test_cvodes_sparse():

    times = np.array([40.0, 400.0])
    x0 = np.array([1.0, 0.0, 0.0])
    nnz = 7
    p = np.array([0.04, 1.0e4, 3.0e7])
    k = np.array([1.0])
    options = {"AbsTol": 1e-12, "RelTol": 1e-08, "linsol": "superLU"}

    t, x = suneigen4py.cvodes(ode, jac, x0, times, nnz, p, k, options)

    np.testing.assert_almost_equal(x[0,0], 7.158017e-01, decimal=4)
    np.testing.assert_almost_equal(x[0,1], 9.1850e-06, decimal=9)
    np.testing.assert_almost_equal(x[0,2], 2.8418e-01, decimal=4)
    np.testing.assert_almost_equal(x[1,0], 4.5053e-01, decimal=4)
    np.testing.assert_almost_equal(x[1,1], 3.2232e-06, decimal=9)
    np.testing.assert_almost_equal(x[1,2], 5.4946e-01, decimal=4)
