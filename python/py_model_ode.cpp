#include "py_model_ode.h"

namespace suneigen::pymodel{

    void PyModel_ODE::fxdot(realtype *xdot, const realtype t, const realtype *x) {

        pybind11::array_t<double> x_ = pybind11::array_t<double>(nx, x);

        pybind11::array_t<double> xdot_ = fxdot_cb_(t, x_);

        // Copy result to xdot
        pybind11::buffer_info xdot_buf = xdot_.request();
        auto xdot_ptr = static_cast<double*>(xdot_buf.ptr);
        for (int i = 0; i < nx; ++i) {
            xdot[i] = xdot_ptr[i];
        }
    }

    void PyModel_ODE::fJSparse(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x) {

        pybind11::array_t<double> x_ = pybind11::array_t<double>(nx, x);

        // Evaluate the jacobian
        pybind11::array_t<double, pybind11::array::f_style> fj_ = fJSparse_cb(t, x_);
        pybind11::buffer_info jac_buf = fj_.request();
        auto *jac_ptr = static_cast<double *>(jac_buf.ptr);

        // Copy the jacobian to c++ library (can we swap instead?)
        for (int i = 0; i < nx*nx; ++i) {
            JSparse->data[i] = jac_ptr[i];
        }
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < nx; ++j) {
                JSparse->indexvals[j + i * nx] = j;
            }
        }
        for (int i = 0; i < nx + 1; i++){
            JSparse->indexptrs[i] = i*nx;
        }

        /*
        // Develop sparse format
        Eigen::Map<Eigen::MatrixXd> jac_map(jac_ptr, 3, 3);
        Eigen::SparseMatrix<double> jac_sp = jac_map.sparseView();

        int counter = 0;
        for (int k=0; k < jac_sp.outerSize(); ++k) {
            JSparse->indexptrs[k] = jac_sp.outerIndexPtr()[k];
            for (Eigen::SparseMatrix<double>::InnerIterator it(jac_sp, k); it; ++it){
                JSparse->indexvals[counter] = it.row();
                JSparse->data[counter] = it.value();
                ++counter;
            }
        }
        */
    }

    void PyModel_ODE::fx0(realtype *x0, const realtype t) {
        for (int i = 0; i < nx; ++i) {
            pybind11::buffer_info x0_buf = x0_.request();
            auto *x0_ptr = static_cast<double *>(x0_buf.ptr);
            x0[i] = x0_ptr[i];
        }
    }


} // namespace suneigen