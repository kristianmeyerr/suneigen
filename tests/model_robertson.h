#ifndef SUNEIGEN_MODEL_ROBERTSON_H
#define SUNEIGEN_MODEL_ROBERTSON_H

#include "defines.h"
#include "model_ode.h"
#include "model.h"

namespace suneigen::model_robertson{

    extern void xdot_model_robertson(realtype *xdot, realtype t, const realtype *x);
    extern void JSparse_model_robertson(SUNMatrixContent_Sparse JSparse, realtype t, const realtype *x);
    extern void fx0_model_robertson(realtype* x0, realtype t);

    class Model_robertson : public suneigen::Model_ODE {
    public:

        Model_robertson() : suneigen::Model_ODE(
                suneigen::ModelDimensions(3, 0, 9))
                {}

        void fxdot(realtype *xdot, const realtype t, const realtype *x) override {
            xdot_model_robertson(xdot, t, x);
        }

        void fJSparse(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x) override {
            JSparse_model_robertson(JSparse, t, x);
        }

        void fx0(realtype *x0, const realtype t) override {
            fx0_model_robertson(x0, t);
        }
    };
}

namespace suneigen::generic_model {

        std::unique_ptr<suneigen::Model> getModel();

    } // namespace suneigen

#endif //SUNEIGEN_MODEL_ROBERTSON_H
