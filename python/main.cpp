#include <pybind11/pybind11.h>
#include "pycvodes.h"
#include <pybind11/eigen.h>
#include <Eigen/SparseCore>
#include <pybind11/functional.h>
#include <iostream>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
namespace py = pybind11;


pybind11::array_t<double> jac(const std::function<pybind11::array_t<double>(double, pybind11::array_t<double>)> &f) {
    auto x = py::array_t<double>(3);
    py::buffer_info x_buf = x.request();
    auto *x_ptr = static_cast<double *>(x_buf.ptr);
    x_ptr[0] = 1.0;
    x_ptr[1] = 2.0;
    x_ptr[2] = 3.0;

    pybind11::array_t<double, py::array::f_style> jac = f(0.0, x);
    py::buffer_info jac_buf = jac.request();
    auto *jac_ptr = static_cast<double *>(jac_buf.ptr);

    Eigen::Map<Eigen::MatrixXd> jac_map(jac_ptr, 3, 3);

    std::cout << jac_map << std::endl;

    Eigen::SparseMatrix<double> jac_sp = jac_map.sparseView();

    for (int i = 0; i < 4; ++i) {
        std::cout << jac_sp.outerIndexPtr()[i] << std::endl;
    }

    return jac;

}

PYBIND11_MODULE(suneigen4py, m) {

    m.def("cvodes", &cvodes);

    m.def("jac", &jac);


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
