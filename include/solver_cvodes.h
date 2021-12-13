//
// Created by Kristian Meyer on 04/11/2021.
//

#ifndef SUNEIGEN_SOLVER_CVODES_H
#define SUNEIGEN_SOLVER_CVODES_H


#include "solver.h"
#include "vector.h"

namespace suneigen {

    /**
     * @brief The CVodeSolver class is a wrapper around the SUNDIALS CVODES solver.
     */
    class CVodeSolver : public Solver {
    public:
        using Solver::Solver;

        /**
         * @brief Clone this instance
         * @return The clone
         */
        Solver *clone() const override;

        void reInit(realtype t0, const Vector &yy0,
                    const Vector &yp0) const override;

        void sensReInit(const VectorArray &yyS0,
                        const VectorArray &ypS0) const override;

        void sensToggleOff() const override;

        void getSensDky(realtype t, int k) const override;

        void setStopTime(realtype tstop) const override;

        void turnOffRootFinding() const override;

        int solve(realtype tout, int itask) const override;

        void getRootInfo(int *rootsfound) const override;

        const Model *getModel() const override;

        void setLinearSolver() const override;

        void setNonLinearSolver() const override;

        void setNonLinearSolverSens() const override;

    protected:

        void calcIC(realtype tout1) const override;

        void allocateSolver() const override;

        void setSStolerances(double rtol, double atol) const override;

        void setSensSStolerances(double rtol,
                                 const double *atol) const override;

        void setSensErrCon(bool error_corr) const override;

        void setErrHandlerFn() const override;

        void reInitPostProcessF(realtype tnext) const override;

        /**
         * @brief Postprocessing of the solver memory after a discontinuity
         * @param cv_mem pointer to CVODES solver memory object
         * @param t pointer to integration time
         * @param yout  new state vector
         * @param tout  anticipated next integration timepoint.
         */
        void reInitPostProcess(void *cv_mem, realtype *t, Vector *yout,
                               realtype tout) const;

        void setUserData() const override;

        void setMaxNumSteps(size_t mxsteps) const override;

        void setStabLimDet(int stldet) const override;

        void setId(const Model *model) const override;

        void setSuppressAlg(bool flag) const override;

        /**
         * @brief resetState reset the CVODES solver to restart integration after a rhs discontinuity.
         * @param cv_mem pointer to CVODES solver memory object
         * @param y0 new state vector
         */
        void resetState(void *cv_mem, const_N_Vector y0) const;

        void setSensParams(const realtype *p, const realtype *pbar,
                           const int *plist) const override;

        void getDky(realtype t, int k) const override;

        void getNumSteps(const void* sun_mem, size_t *numsteps) const override;

        void getNumNonlinSolvIters(const void* sun_mem, size_t *numnonlinsolviters) const override;

        void getNumJacEvals(const void* sun_mem, size_t *numjacevals) const override;

        void getNumRhsEvals(const void* sun_mem,
                            size_t *numrhsevals) const override;

        void getNumErrTestFails(const void* sun_mem,
                                size_t *numerrtestfails) const override;

        void getNumNonlinSolvConvFails(const void* sun_mem,
                                       size_t *numnonlinsolvconvfails) const override;

        void getLastOrder(const void* sun_mem, size_t *order) const override;

        std::string getSolverType() const override;

        void init(realtype t0, const Vector &x0,
                  const Vector &dx0) const override;

        void rootInit(int ne) const override;

        void setDenseJacFn() const override;

        void setSparseJacFn() const override;

        void getSens() const override;

        void getQuad(realtype &t) const override;

        void getQuadDky(realtype t, int k) const override;



    };

}

#endif //SUNEIGEN_SOLVER_CVODES_H
