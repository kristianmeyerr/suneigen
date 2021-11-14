#ifndef SUNEIGEN_PY_MODEL_ODE_H
#define SUNEIGEN_PY_MODEL_ODE_H

#include "defines.h"
#include "model_ode.h"

#include <utility>
#include "model.h"

namespace suneigen::pymodel{

    class PyModel_ODE : public suneigen::Model_ODE {
    public:

        PyModel_ODE(std::function<double(double, double)>  fxdot)
        : suneigen::Model_ODE(
                suneigen::ModelDimensions(3, 0, 9)),
                fxdot_(std::move(fxdot))
                {}

        void fxdot(realtype *xdot, const realtype t, const realtype *x) override {

        }

        void fJSparse(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x) override {
            (void)JSparse;
            (void)t;
            (void)x;
        }

        void fx0(realtype *x0, const realtype t) override {
            (void)t;
            (void)x0;
        }

    protected:

        const std::function<double(double, double)> fxdot_;

    };
}

#endif //SUNEIGEN_PY_MODEL_ODE_H
