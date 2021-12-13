#include "model_ode.h"
#include "solver_cvodes.h"

namespace suneigen {

    void Model_ODE::fJ(const realtype t, const realtype cj, const Vector &x,
                       const Vector & dx, const Vector &xdot,
                       SUNMatrix J) {
        (void)dx; (void) cj;
        fJ(t, x.getNVector(), xdot.getNVector(), J);
    }

    void Model_ODE::fJ(realtype t, const_N_Vector x, const_N_Vector xdot, SUNMatrix J) {
        (void)xdot;
        fJSparse(t, x, derived_state_.J_.get());
        derived_state_.J_.refresh();
        auto JDense = SUNMatrixWrapper(J);
        derived_state_.J_.to_dense(JDense);
    }

    void Model_ODE::fxdot(const realtype t, const Vector &x,
                          const Vector& dx, Vector &xdot) {
        (void)dx;
        fxdot(t, x.getNVector(), xdot.getNVector());
    }

    void Model_ODE::fxdot(realtype t, const_N_Vector x, N_Vector xdot) {
        (void) t;
        (void) x;
        (void) xdot;
        auto x_pos = computeX_pos(x);
        N_VConst(0.0, xdot);
        fxdot(N_VGetArrayPointer(xdot), t,
              N_VGetArrayPointerConst(x_pos),state_.parameters.data(),
              state_.fixedParameters.data(),
              state_.h.data());

    }

    void Model_ODE::froot(const realtype t, const Vector &x,
                          const Vector & /*dx*/, gsl::span<realtype> root) {
        froot(t, x.getNVector(), root);
    }

    void Model_ODE::froot(realtype t, const_N_Vector x, gsl::span<realtype> root) {
        auto x_pos = computeX_pos(x);
        std::fill(root.begin(), root.end(), 0.0);
        froot(root.data(), t, N_VGetArrayPointerConst(x_pos),
              state_.parameters.data(), state_.fixedParameters.data(), state_.h.data());
    }

    void Model_ODE::froot(realtype * /*root*/, const realtype /*t*/,
                          const realtype * /*x*/, const realtype * /*p*/,
                          const realtype * /*k*/, const realtype * /*h*/) {
        throw SunException("Requested functionality is not supported as %s is not "
                           "implemented for this model!",
                           __func__); // not implemented
    }

    void Model_ODE::fJSparse(const realtype t, const realtype /*cj*/,
                             const Vector &x, const Vector & /*dx*/,
                             const Vector & /*xdot*/, SUNMatrix J) {
        fJSparse(t, x.getNVector(), J);
    }

    std::unique_ptr<Solver> Model_ODE::getSolver() {
        return std::unique_ptr<Solver>(new suneigen::CVodeSolver());
    }

    void Model_ODE::fJSparse(realtype t, const_N_Vector x, SUNMatrix J) {

        auto x_pos = computeX_pos(x);
        fJSparse(static_cast<SUNMatrixContent_Sparse>(SM_CONTENT_S(J)), t,
                 N_VGetArrayPointerConst(x_pos), state_.parameters.data(),
                 state_.fixedParameters.data(), state_.h.data());

    }

    void Model_ODE::fJSparse(SUNMatrixContent_Sparse /*JSparse*/,
                             const realtype /*t*/, const realtype * /*x*/,
                             const realtype * /*p*/, const realtype * /*k*/,
                             const realtype * /*h*/) {
        throw SunException("Requested functionality is not supported as %s "
                           "is not implemented for this model!",
                           __func__); // not implemented
    }

}
