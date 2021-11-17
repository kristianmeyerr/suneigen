#ifndef SUNEIGEN_PY_MODEL_ODE_H
#define SUNEIGEN_PY_MODEL_ODE_H

#include "defines.h"
#include "model_ode.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <utility>
#include "model.h"
#include <iostream>
#include <Eigen/SparseCore>

namespace suneigen::pymodel{

    class PyModel_ODE : public Model_ODE {
    public:

        PyModel_ODE(
                std::function<pybind11::array_t<double>(double, pybind11::array_t<double>)> fxdot_cb,
                std::function<pybind11::array_t<double>(double, pybind11::array_t<double>)> fJSparse_cb,
                pybind11::array_t<double> x0)
                : suneigen::Model_ODE(suneigen::ModelDimensions(3, 0, 9)),
                fxdot_cb_(std::move(fxdot_cb)),
                fJSparse_cb(std::move(fJSparse_cb)),
                x0_(std::move(x0))
        {}

        void fxdot(realtype *xdot, realtype t, const realtype *x) override;

        void fJSparse(SUNMatrixContent_Sparse JSparse, realtype t, const realtype *x) override;

        void fx0(realtype *x0, realtype t) override;

    private:
        pybind11::array_t<double> x0_;
        std::function<pybind11::array_t<double>(double, pybind11::array_t<double>)> fxdot_cb_;
        std::function<pybind11::array_t<double>(double, pybind11::array_t<double>)> fJSparse_cb;

    };
}

#endif //SUNEIGEN_PY_MODEL_ODE_H
