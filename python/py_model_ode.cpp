#include "py_model_ode.h"
#include <iostream>
using namespace pybind11::literals;

namespace suneigen::pymodel{

    void PyModel_ODE::fxdot(realtype *xdot, const realtype t, const realtype *x,
                            const realtype *p, const realtype *k, const realtype *h) {

        pybind11::array_t<double> x_ = pybind11::array_t<double>(static_cast<ssize_t>(nx), x);
        pybind11::array_t<double> p_ = pybind11::array_t<double>(static_cast<ssize_t>(np()), x);
        pybind11::array_t<double> k_ = pybind11::array_t<double>(static_cast<ssize_t>(nk()), x);
        pybind11::array_t<double> h_ = pybind11::array_t<double>(ne, h);

        pybind11::array_t<double> xdot_ = fxdot_cb_(t, x_, p_, k_, h_);

        // Copy result to xdot
        pybind11::buffer_info xdot_buf = xdot_.request();
        auto xdot_ptr = static_cast<double*>(xdot_buf.ptr);
        for (int i = 0; i < nx; ++i) {
            xdot[i] = xdot_ptr[i];
        }
    }

    void PyModel_ODE::froot(realtype *root, const realtype t, const realtype *x,
               const realtype *p, const realtype *k, const realtype *h){
    }

    void PyModel_ODE::fJSparse(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x,
                               const realtype *p, const realtype *k, const realtype *h) {

        pybind11::array_t<double> x_ = pybind11::array_t<double>(static_cast<ssize_t>(nx), x);
        pybind11::array_t<double> p_ = pybind11::array_t<double>(static_cast<ssize_t>(np()), x);
        pybind11::array_t<double> k_ = pybind11::array_t<double>(static_cast<ssize_t>(nk()), x);
        pybind11::array_t<double> h_ = pybind11::array_t<double>(ne, x);

        // Evaluate the jacobian
        pybind11::tuple result = fJSparse_cb(t, x_, p_, k_, h_);

        // Iterate through the tuple results
        auto result_ptr = result.begin();

        auto indexvals = pybind11::cast<pybind11::array_t<sunindextype>>(*result_ptr);
        pybind11::buffer_info indexvals_buf = indexvals.request();
        auto indexvals_ptr = static_cast<sunindextype*>(indexvals_buf.ptr);
        for (int i = 0; i < indexvals.size(); ++i) {
            JSparse->indexvals[i] = indexvals_ptr[i];
        }

        result_ptr++;
        auto indexptrs = pybind11::cast<pybind11::array_t<sunindextype>>(*result_ptr);
        pybind11::buffer_info indexptrs_buf = indexptrs.request();
        auto indexptrs_ptr = static_cast<sunindextype*>(indexptrs_buf.ptr);
        for (int i = 0; i < indexptrs.size(); ++i) {
            JSparse->indexptrs[i] = indexptrs_ptr[i];
        }

        result_ptr++;
        auto data = pybind11::cast<pybind11::array_t<double>>(*result_ptr);
        pybind11::buffer_info data_buf = data.request();
        auto data_ptr = static_cast<double*>(data_buf.ptr);
        for (int i = 0; i < data.size(); ++i) {
            JSparse->data[i] = data_ptr[i];
        }

    }

    void PyModel_ODE::fx0(realtype *x0, const realtype /*t*/, const realtype */*p*/, const realtype */*k*/) {
        for (int i = 0; i < nx; ++i) {
            pybind11::buffer_info x0_buf = x0_.request();
            auto *x0_ptr = static_cast<double *>(x0_buf.ptr);
            x0[i] = x0_ptr[i];
        }
    }


} // namespace suneigen