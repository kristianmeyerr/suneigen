#ifndef SUNEIGEN_PYCVODES_H
#define SUNEIGEN_PYCVODES_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include "py_model_ode.h"

namespace suneigen {

    pybind11::tuple cvodes(const pyfunction_cb &ode,
                           const pyjac_cb &jac,
                           const pybind11::array_t<double>& x0,
                           const pybind11::array_t<double>& times,
                           const int& nnz,
                           const pybind11::array_t<double>& p,
                           const pybind11::array_t<double>& k,
                           const pybind11::dict& options);

}

#endif //SUNEIGEN_PYCVODES_H
