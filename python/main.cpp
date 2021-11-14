#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include "suneigen.h"


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
/*
py::tuple cvodes(const std::function<py::array_t<double>(double, py::array_t<double>)> &ode,
                 const std::function<py::array_t<double>(double)> &ode_x0, const py::array_t<double>& times) {

    py::buffer_info buf1 = times.request();
    auto times_ptr = static_cast<double *>(buf1.ptr);

    auto time = py::array_t<double>(3);
    py::buffer_info time_buffer = time.request();
    auto *ptr1 = static_cast<realtype *>(time_buffer.ptr);
    ptr1[0] = times_ptr[0];
    ptr1[1] = times_ptr[1];
    ptr1[2] = times_ptr[2];

    py::array_t<double> x0 = ode_x0(0.0);

    return py::make_tuple(time, ode(0.0, x0));

}
*/

int add(int a, int b);
int add(int a, int b){
    return a + b + 2;
}

PYBIND11_MODULE(suneigen4py, m) {
    m.def("add", &add);


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
