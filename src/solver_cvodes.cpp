#include "solver_cvodes.h"

#include <model_ode.h>
#include "sundials_matrix_wrapper.h"

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_diag.h>
#include <cvodes/cvodes_impl.h>

#define ZERO RCONST(0.0)
#define ONE RCONST(1.0)

namespace suneigen {

    /*
     * The following static members are callback function to CVODES.
     * Their signatures must not be changes.
     */
    static int fxdot(realtype t, N_Vector x, N_Vector xdot, void* user_data);

    static int fJ(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                  void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    static int fJSparse(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                        void *user_data, N_Vector tmp1, N_Vector tmp2,
                        N_Vector tmp3);

    static int froot(realtype t, N_Vector x, realtype *root, void *user_data);

    // Ensure SunEigen options are in sync with Sundials options
    static_assert(static_cast<int>(InternalSensitivityMethod::simultaneous) == CV_SIMULTANEOUS);

    Solver* CVodeSolver::clone() const { return new CVodeSolver(*this); }

    void CVodeSolver::sensToggleOff() const {
        auto status = CVodeSensToggleOff(solver_memory_.get());
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSensToggleOff");
    }

    void CVodeSolver::reInit(const realtype t0, const Vector &yy0,
                             const Vector & /*yp0*/) const {
        auto cv_mem = static_cast<CVodeMem>(solver_memory_.get());
        cv_mem->cv_tn = t0;
        if (solver_was_called_F_)
            force_reinit_postprocess_F_ = true;
        x_.copy(yy0);
        resetState(cv_mem, x_.getNVector());
    }

    void CVodeSolver::sensReInit(const VectorArray &yyS0,
                                 const VectorArray & /*ypS0*/) const {
        auto cv_mem = static_cast<CVodeMem>(solver_memory_.get());
        /* Initialize znS[0] in the history array */
        for (unsigned int is = 0; is < np(); is++)
            cv_mem->cv_cvals[is] = ONE;
        if (solver_was_called_F_)
            force_reinit_postprocess_F_ = true;
        sx_.copy(yyS0);
        int status = N_VScaleVectorArray(static_cast<int>(np()), cv_mem->cv_cvals,
                                         sx_.getNVectorArray(), cv_mem->cv_znS[0]);
        if (status != CV_SUCCESS)
            throw CvodeException(CV_VECTOROP_ERR, "CVodeSensReInit");
    }

    void CVodeSolver::allocateSolver() const {
        if (!solver_memory_)
            solver_memory_ = std::unique_ptr<void, std::function<void(void *)>>(
                    CVodeCreate(static_cast<int>(lmm_)),
                    [](void *ptr) { CVodeFree(&ptr); });
    }

    void CVodeSolver::getNumSteps(const void *sun_mem, size_t* numsteps) const {
        int status = CVodeGetNumSteps(const_cast<void *>(sun_mem), reinterpret_cast<long *>(numsteps));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetNumSteps");
    }

    void CVodeSolver::getNumRhsEvals(const void *sun_mem,
                                     size_t *numrhsevals) const {
        int status = CVodeGetNumRhsEvals(const_cast<void *>(sun_mem), reinterpret_cast<long *>(numrhsevals));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetNumRhsEvals");
    }

    void CVodeSolver::getNumErrTestFails(const void *sun_mem,
                                         size_t *numerrtestfails) const {
        int status =
                CVodeGetNumErrTestFails(const_cast<void *>(sun_mem), reinterpret_cast<long *>(numerrtestfails));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetNumErrTestFails");
    }

    void CVodeSolver::getNumNonlinSolvConvFails(
            const void *sun_mem, size_t *numnonlinsolvconvfails) const {
        int status = CVodeGetNumNonlinSolvConvFails(const_cast<void *>(sun_mem),
                                                    reinterpret_cast<long *>(numnonlinsolvconvfails));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetNumNonlinSolvConvFails");
    }

    void CVodeSolver::getNumNonlinSolvIters(const void *sun_mem, size_t *numnonlinsolviters) const {
        int status = CVodeGetNumNonlinSolvIters(const_cast<void *>(sun_mem),
                                                reinterpret_cast<long *>(numnonlinsolviters));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "getNumNonlinSolvIters");
    }

    void CVodeSolver::getNumJacEvals(const void* sun_mem, size_t *numjacevals) const {
        int status = CVodeGetNumJacEvals(const_cast<void *>(sun_mem),
                                                reinterpret_cast<long *>(numjacevals));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "getNumJacEvals");
    }

    std::string CVodeSolver::getSolverType() const {
        return std::string("CVode");
    }

    void CVodeSolver::getLastOrder(const void *sun_mem, size_t *order) const {
        int status = CVodeGetLastOrder(const_cast<void *>(sun_mem), reinterpret_cast<int *>(order));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetLastOrder");
    }

    void CVodeSolver::setStopTime(const realtype tstop) const {
        int status = CVodeSetStopTime(solver_memory_.get(), tstop);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetStopTime");
    }

    void CVodeSolver::turnOffRootFinding() const {
        int status = CVodeRootInit(solver_memory_.get(), 0, nullptr);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeRootInit");
    }

    void CVodeSolver::init(const realtype t0, const Vector &x0,
                           const Vector & dx0) const {

        (void) dx0;  // Not needed for cvodes
        solver_was_called_F_ = false;
        force_reinit_postprocess_F_ = false;
        t_ = t0;
        x_ = x0;
        int status;
        if (getInitDone()) {
            status = CVodeReInit(solver_memory_.get(), t0, x_.getNVector());
        } else {
            status = CVodeInit(solver_memory_.get(), fxdot, t0, x_.getNVector());
            setInitDone();
        }
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeInit");
    }

    void CVodeSolver::calcIC(const realtype /*tout1*/) const {}

    void CVodeSolver::rootInit(int ne) const {
        int status = CVodeRootInit(solver_memory_.get(), ne, froot);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeRootInit");
    }

    void CVodeSolver::setDenseJacFn() const {
        int status = CVodeSetJacFn(solver_memory_.get(), fJ);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetJacFn");
    }

    void CVodeSolver::setUserData() const {
        int status = CVodeSetUserData(solver_memory_.get(), &user_data);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetUserData");
    }

    void CVodeSolver::setMaxNumSteps(size_t mxsteps) const {
        int status = CVodeSetMaxNumSteps(solver_memory_.get(), static_cast<long>(mxsteps));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetMaxNumSteps");
    }

    void CVodeSolver::reInitPostProcessF(const realtype tnext) const {
        reInitPostProcess(solver_memory_.get(), &t_, &x_, tnext);
        force_reinit_postprocess_F_ = false;
    }

    void CVodeSolver::reInitPostProcess(void *ami_mem, realtype *t, Vector *yout,
                                        const realtype tout) const {
        auto cv_mem = static_cast<CVodeMem>(ami_mem);
        auto nst_tmp = cv_mem->cv_nst;
        cv_mem->cv_nst = 0;

        auto status = CVodeSetStopTime(cv_mem, tout);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetStopTime");

        status = CVode(ami_mem, tout, yout->getNVector(), t, CV_ONE_STEP);

        if (status == CV_ROOT_RETURN)
            throw CvodeException(status, "CVode returned a root after "
                                         "reinitialization. The initial step-size after the event or "
                                         "heaviside function is too small. To fix this, increase absolute "
                                         "and relative tolerances!");
        if (status != CV_SUCCESS)
            throw CvodeException(status, "reInitPostProcess");

        cv_mem->cv_nst = nst_tmp + 1;
        if (cv_mem->cv_adjMallocDone == SUNTRUE) {
            /* add new step to history array, this is copied from CVodeF */
            auto ca_mem = cv_mem->cv_adj_mem;
            auto dt_mem = ca_mem->dt_mem;

            if (cv_mem->cv_nst % ca_mem->ca_nsteps == 0) {
                /* currently not implemented, we should never get here as we
                 limit cv_mem->cv_nst < ca_mem->ca_nsteps, keeping this for
                 future regression */
                throw CvodeException(SUNEIGEN_ERROR, "reInitPostProcess");
            }

            /* Load next point in dt_mem */
            dt_mem[cv_mem->cv_nst % ca_mem->ca_nsteps]->t = *t;
            ca_mem->ca_IMstore(cv_mem, dt_mem[cv_mem->cv_nst % ca_mem->ca_nsteps]);

            /* Set t1 field of the current ckeck point structure
             for the case in which there will be no future
             check points */
            ca_mem->ck_mem->ck_t1 = *t;

            /* tfinal is now set to *tret */
            ca_mem->ca_tfinal = *t;
        }
    }

    int CVodeSolver::solve(const realtype tout, const int itask) const {
        if (force_reinit_postprocess_F_)
            reInitPostProcessF(tout);
        int status = CVode(solver_memory_.get(), tout, x_.getNVector(), &t_, itask);
        if (status < 0) // status > 0 is okay and is used for e.g. root return
            throw IntegrationFailure(status, t_);
        solver_was_called_F_ = true;
        return status;
    }

    void CVodeSolver::setStabLimDet(const int stldet) const {
        int status = CVodeSetStabLimDet(solver_memory_.get(), stldet);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetStabLimDet");
    }

    void CVodeSolver::setId(const Model * /*model*/) const {}

    void CVodeSolver::setSuppressAlg(const bool /*flag*/) const {}

    void CVodeSolver::resetState(void *ami_mem, const_N_Vector y0) const {

        auto cv_mem = static_cast<CVodeMem>(ami_mem);
        /* here we force the order in the next step to zero, and update the
         Nordsieck history array, this is largely copied from CVodeReInit with
         explanations from cvodes_impl.h
         */

        /* Set step parameters */

        /* current order */
        cv_mem->cv_q = 1;
        /* L = q + 1 */
        cv_mem->cv_L = 2;
        /* number of steps to wait before updating in q */
        cv_mem->cv_qwait = cv_mem->cv_L;
        /* last successful q value used                */
        cv_mem->cv_qu = 0;
        /* last successful h value used                */
        cv_mem->cv_hu = ZERO;
        /* tolerance scale factor                      */
        cv_mem->cv_tolsf = ONE;

        /* Initialize other integrator optional outputs */

        /* actual initial stepsize                     */
        cv_mem->cv_h0u = ZERO;
        /* step size to be used on the next step       */
        cv_mem->cv_next_h = ZERO;
        /* order to be used on the next step           */
        cv_mem->cv_next_q = 0;

        /* write updated state to Nordsieck history array  */
        N_VScale(ONE, const_cast<N_Vector>(y0), cv_mem->cv_zn[0]);
    }

    void CVodeSolver::setSensParams(const realtype *p, const realtype *pbar,
                                    const int *plist) const {
        int status = CVodeSetSensParams(
                solver_memory_.get(), const_cast<realtype *>(p),
                const_cast<realtype *>(pbar), const_cast<int *>(plist));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetSensParams");
    }

    void CVodeSolver::setSparseJacFn() const {
        int status = CVodeSetJacFn(solver_memory_.get(), fJSparse);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetJacFn");
    }

    void CVodeSolver::getRootInfo(int *rootsfound) const {
        int status = CVodeGetRootInfo(solver_memory_.get(), rootsfound);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetRootInfo");
    }

    void CVodeSolver::setLinearSolver() const {
        int status = CVodeSetLinearSolver(solver_memory_.get(), linear_solver_->get(),
                                          linear_solver_->getMatrix());
        if (status != CV_SUCCESS)
            throw CvodeException(status, "setLinearSolver");
    }

    void CVodeSolver::setNonLinearSolver() const {
        int status =
                CVodeSetNonlinearSolver(solver_memory_.get(), non_linear_solver_->get());
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetNonlinearSolver");
    }

    void CVodeSolver::setNonLinearSolverSens() const {
        if (getSensitivityOrder() < SensitivityOrder::first)
            return;
        if (getSensitivityMethod() != SensitivityMethod::forward)
            return;

        int status = CV_SUCCESS;

        assert(ism_==InternalSensitivityMethod::staggered || ism_==InternalSensitivityMethod::simultaneous
        || ism_ ==InternalSensitivityMethod::staggered1);

        switch (ism_) {
            case InternalSensitivityMethod::staggered:
                status = CVodeSetNonlinearSolverSensStg(solver_memory_.get(),
                                                        non_linear_solver_sens_->get());
                break;
            case InternalSensitivityMethod::simultaneous:
                status = CVodeSetNonlinearSolverSensSim(solver_memory_.get(),
                                                        non_linear_solver_sens_->get());
                break;
            case InternalSensitivityMethod::staggered1:
                status = CVodeSetNonlinearSolverSensStg1(solver_memory_.get(),
                                                         non_linear_solver_sens_->get());
                break;
        }

        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSolver::setNonLinearSolverSens");
    }

    void CVodeSolver::setSStolerances(const double rtol, const double atol) const {
        int status = CVodeSStolerances(solver_memory_.get(), rtol, atol);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSStolerances");
    }

    void CVodeSolver::setSensSStolerances(const double rtol,
                                          const double *atol) const {
        int status = CVodeSensSStolerances(solver_memory_.get(), rtol,
                                           const_cast<double *>(atol));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSensEEtolerances");
    }

    void CVodeSolver::setErrHandlerFn() const {
        int status =
                CVodeSetErrHandlerFn(solver_memory_.get(), wrapErrHandlerFn,
                                     reinterpret_cast<void*>(
                                             const_cast<CVodeSolver*>(this)));
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetErrHandlerFn");
    }

    void CVodeSolver::getDky(realtype t, int k) const {
        int status = CVodeGetDky(solver_memory_.get(), t, k, dky_.getNVector());
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetDky");
    }

    void CVodeSolver::getSens() const {
        realtype tDummy = 0;
        int status =
            CVodeGetSens(solver_memory_.get(), &tDummy, sx_.getNVectorArray());
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetSens");
    }

    void CVodeSolver::getQuad(realtype &t) const {
        int status = CVodeGetQuad(solver_memory_.get(), &t, xQ_.getNVector());
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetQuad");
    }

    void CVodeSolver::getQuadDky(const realtype t, const int k) const {
        int status = CVodeGetQuadDky(solver_memory_.get(), t, k, xQ_.getNVector());
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetQuadDky");
    }

    void CVodeSolver::getSensDky(const realtype t, const int k) const {
        int status =
                CVodeGetSensDky(solver_memory_.get(), t, k, sx_.getNVectorArray());
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeGetSens");
    }

    void CVodeSolver::setSensErrCon(const bool error_corr) const {
        int status = CVodeSetSensErrCon(solver_memory_.get(), error_corr);
        if (status != CV_SUCCESS)
            throw CvodeException(status, "CVodeSetSensErrCon");
    }

    /**
     * @brief Event trigger function for events
     * @param t timepoint
     * @param x Vector with the states
     * @param root array with root function values
     * @param user_data object with user input
     * @return status flag indicating successful execution
     */
    static int froot(realtype t, N_Vector x, realtype *root,
                     void *user_data) {
        auto typed_udata = static_cast<CVodeSolver::user_data_type *>(user_data);
        Expects(typed_udata);
        auto model = dynamic_cast<Model_ODE *>(typed_udata->first);
        Expects(model);

        model->froot(t, x, gsl::make_span<realtype>(root, model->ne));
        return model->checkFinite(gsl::make_span<realtype>(root, model->ne),
                                  "root function");
    }

    /**
     * @brief residual function of the ODE
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     * @param user_data object with user input
     * @return status flag indicating successful execution
     */
    static int fxdot(realtype t, N_Vector x, N_Vector xdot, void *user_data) {
        auto typed_udata = static_cast<CVodeSolver::user_data_type *>(user_data);
        Expects(typed_udata);
        auto model = dynamic_cast<Model_ODE *>(typed_udata->first);
        Expects(model);
        auto solver = dynamic_cast<CVodeSolver const *>(typed_udata->second);
        Expects(model);

        if(solver->timeExceeded()) {
            return SUNEIGEN_MAX_TIME_EXCEEDED;
        }

        if (t > 1e200 && !model->checkFinite(gsl::make_span(x), "fxdot")) {
            /* when t is large (typically ~1e300), CVODES may pass all NaN x
               to fxdot from which we typically cannot recover. To save time
               on normal execution, we do not always want to check finiteness
               of x, but only do so when t is large and we expect problems. */
            return SUNEIGEN_UNRECOVERABLE_ERROR;
        }

        model->fxdot(t, x, xdot);
        return model->checkFinite(gsl::make_span(xdot), "fxdot");
    }

    /**
     * @brief Jacobian of xdot with respect to states x
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written
     * @param user_data object with user input
     * @param tmp1 temporary storage vector
     * @param tmp2 temporary storage vector
     * @param tmp3 temporary storage vector
     * @return status flag indicating successful execution
     **/
    static int fJ(realtype t, N_Vector x, N_Vector xdot, SUNMatrix J,
                  void *user_data, N_Vector tmp1, N_Vector tmp2,
                  N_Vector tmp3) {
        (void)tmp1; (void)tmp2; (void)tmp3;

        auto typed_udata = static_cast<CVodeSolver::user_data_type *>(user_data);
        Expects(typed_udata);
        auto model = dynamic_cast<Model_ODE *>(typed_udata->first);
        Expects(model);

        model->fJ(t, x, xdot, J);
        return model->checkFinite(gsl::make_span(J), "Jacobian");
    }

    /**
     * @brief J in sparse form (for sparse solvers from the SuiteSparse Package)
     * @param t timepoint
     * @param x Vector with the states
     * @param xdot Vector with the right hand side
     * @param J Matrix to which the Jacobian will be written
     * @param user_data object with user input
     * @param tmp1 temporary storage vector
     * @param tmp2 temporary storage vector
     * @param tmp3 temporary storage vector
     * @return status flag indicating successful execution
     */
    static int fJSparse(realtype t, N_Vector x, N_Vector xdot,
                        SUNMatrix J, void *user_data, N_Vector tmp1,
                        N_Vector tmp2, N_Vector tmp3) {
        (void)tmp1; (void)tmp2; (void)tmp3; (void)xdot;

        auto typed_udata = static_cast<CVodeSolver::user_data_type *>(user_data);
        Expects(typed_udata);
        auto model = dynamic_cast<Model_ODE *>(typed_udata->first);
        Expects(model);

        model->fJSparse(t, x, J);
        return model->checkFinite(gsl::make_span(J), "Jacobian");
    }

    const Model *CVodeSolver::getModel() const {
        if (!solver_memory_)
            throw SunException("Solver has not been allocated, information is not "
                               "available");
        auto cv_mem = static_cast<CVodeMem>(solver_memory_.get());

        auto typed_udata = static_cast<user_data_type *>(cv_mem->cv_user_data);
        Expects(typed_udata);
        return typed_udata->first;
    }



}
