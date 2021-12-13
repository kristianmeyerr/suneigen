#ifndef SUNEIGEN_MODEL_ROBERTSON_H
#define SUNEIGEN_MODEL_ROBERTSON_H

#include "defines.h"
#include "model_ode.h"
#include "model.h"

#include <vector>

namespace suneigen::model_robertson{

    extern void xdot_model_robertson(realtype *xdot, realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    extern void JSparse_model_robertson(SUNMatrixContent_Sparse JSparse, realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    extern void fx0_model_robertson(realtype* x0, realtype t, const realtype *p, const realtype *k);
    extern void root_model_robertson(realtype *root, realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);
    extern void fz_model_robertson(double *z, unsigned int ie, realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h);


    class Model_robertson : public suneigen::Model_ODE {
    public:

        Model_robertson() : suneigen::Model_ODE(
                suneigen::ModelDimensions(
                        3,
                        3,
                        1,
                        0,
                        2,
                        2,
                        9,
                        0,
                        0),
                suneigen::SimulationParameters(
                        std::vector<realtype>(1, 1.0),
                        std::vector<realtype>(3, 1.0)),
                std::vector<int>{1, 2})
                {}

        void froot(realtype *root, const realtype t, const realtype *x,
                   const realtype *p, const realtype *k, const realtype *h) override {
            root_model_robertson(root, t, x, p, k, h);
        }

        void fxdot(realtype *xdot, const realtype t, const realtype *x,
                const realtype *p, const realtype *k, const realtype *h) override {
            xdot_model_robertson(xdot, t, x, p, k, h);
        }

        void fJSparse(SUNMatrixContent_Sparse JSparse, const realtype t, const realtype *x,
                const realtype *p, const realtype *k, const realtype *h) override {
            JSparse_model_robertson(JSparse, t, x, p, k, h);
        }

        void fx0(realtype *x0, const realtype t, const realtype *p, const realtype *k) override {
            fx0_model_robertson(x0, t, p, k);
        }

        void fsx0(realtype */*sx0*/, const realtype /*t*/,const realtype */*x0*/, const realtype */*p*/, const unsigned int /*ip*/) override {
        }

        void fdeltax(double */*deltax*/, const realtype /*t*/, const realtype */*x*/, const realtype */*p*/,
                const realtype */*h*/, const unsigned int /*ie*/, const realtype */*xdot*/, const realtype */*xdot_old*/) override {
        }

        void fstau(double */*stau*/, const realtype /*t*/, const realtype */*x*/, const realtype */*p*/,
                const realtype */*h*/, const realtype */*sx*/, const unsigned int /*ip*/, const unsigned int /*ie*/) override {
        }

        void fz(double *z, const unsigned int ie, const realtype t, const realtype *x,
                const realtype *p, const realtype *k, const realtype *h) override {
            fz_model_robertson(z, ie, t, x, p, k, h);
        }

    };
}

namespace suneigen::generic_model {

    std::unique_ptr<suneigen::Model> getModel();

} // namespace suneigen

#endif //SUNEIGEN_MODEL_ROBERTSON_H
