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

namespace suneigen {

    using pyfunction_cb = std::function<pybind11::array_t<double>(double, pybind11::array_t<double>,
            pybind11::array_t<double>, pybind11::array_t<double>, pybind11::array_t<double>)>;
    using pyjac_cb = std::function<pybind11::tuple(double, pybind11::array_t<double>,
            pybind11::array_t<double>, pybind11::array_t<double>, pybind11::array_t<double>)>;

    namespace pymodel {

    class PyModel_ODE : public Model_ODE {
        public:

            PyModel_ODE(
                    size_t nx,
                    size_t np,
                    size_t nk,
                    size_t ny,
                    size_t nz,
                    size_t ne,
                    size_t nnz,
                    size_t ubw,
                    size_t lbw,
                    pyfunction_cb fxdot_cb,
                    pyjac_cb fJSparse_cb,
                    pybind11::array_t<double> x0)
                    : suneigen::Model_ODE(
                    suneigen::ModelDimensions(
                            nx,
                            np,
                            nk,
                            ny,
                            nz,
                            ne,
                            nnz,
                            ubw,
                            lbw),
                    suneigen::SimulationParameters(
                            std::vector<realtype>(nk, 1.0),
                            std::vector<realtype>(np, 1.0)),
                    std::vector<int>{}),
                      fxdot_cb_(std::move(fxdot_cb)),
                      fJSparse_cb(std::move(fJSparse_cb)),
                      x0_(std::move(x0)) {}

            void froot(realtype *root, realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h) override;

            void fxdot(realtype *xdot, realtype t, const realtype *x,
                       const realtype *p, const realtype *k, const realtype *h) override;

            void fJSparse(SUNMatrixContent_Sparse JSparse, realtype t, const realtype *x,
                          const realtype *p, const realtype *k, const realtype *h) override;

            void fx0(realtype *x0, realtype t, const realtype *p, const realtype *k) override;

        private:
            pybind11::array_t<double> x0_;
            pyfunction_cb fxdot_cb_;
            pyjac_cb fJSparse_cb;
        };
    }  // namespace pymodel
}  // namespace suneigen

#endif //SUNEIGEN_PY_MODEL_ODE_H
