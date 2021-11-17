#ifndef SUNEIGEN_PYCVODES_H
#define SUNEIGEN_PYCVODES_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

pybind11::tuple cvodes(const std::function<pybind11::array_t<double>(double, pybind11::array_t<double>)> &ode,
                       const std::function<pybind11::array_t<double>(double, pybind11::array_t<double>)> &jac,
                       const pybind11::array_t<double>& x0, const pybind11::array_t<double>& times,
                       const pybind11::dict& options);

#endif //SUNEIGEN_PYCVODES_H
