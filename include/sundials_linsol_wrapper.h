#ifndef SUNEIGEN_SUNDIALS_LINSOL_WRAPPER_H
#define SUNEIGEN_SUNDIALS_LINSOL_WRAPPER_H

#include "exception.h"
#include "sundials_matrix_wrapper.h"
#include "vector.h"

#include <sundials/sundials_config.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol_superlu.h>

#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

namespace suneigen {

    /**
     * @brief A RAII wrapper for SUNLinearSolver structs.
     *
     * For details on member functions see documentation in
     * sunlinsol/sundials_linearsolver.h.
     */
    class SUNLinSolWrapper {
    public:
        SUNLinSolWrapper() = default;

        /**
         * @brief Wrap existing SUNLinearSolver
         * @param linsol The linear solver
         */
        explicit SUNLinSolWrapper(SUNLinearSolver linsol);

        virtual ~SUNLinSolWrapper();

        /**
         * @brief Copy constructor
         * @param other d
         */
        SUNLinSolWrapper(const SUNLinSolWrapper &other) = delete;

        /**
         * @brief Move constructor
         * @param other d
         */
        SUNLinSolWrapper(SUNLinSolWrapper &&other) noexcept;

        /**
         * @brief Copy assignment
         * @param other d
         * @return d
         */
        SUNLinSolWrapper &operator=(const SUNLinSolWrapper &other) = delete;

        /**
         * @brief Move assignment
         * @param other d
         * @return d
         */
        SUNLinSolWrapper &operator=(SUNLinSolWrapper &&other) noexcept = delete;

        /**
         * @brief Returns the wrapped SUNLinSol.
         * @return SUNLinearSolver
         */
        [[nodiscard]] SUNLinearSolver get() const;

        /**
         * @brief Returns an identifier for the linear solver type.
         * @return d
         */
        [[nodiscard]] SUNLinearSolver_Type getType() const;

        /**
         * @brief Performs any linear solver setup needed, based on an updated
         * system matrix A.
         * @param A d
         */
        void setup(SUNMatrix A) const;

        /**
         * @brief Performs any linear solver setup needed, based on an updated
         * system matrix A.
         * @param A d
         */
        void setup(const SUNMatrixWrapper& A) const;

        /**
         * @brief Solves a linear system A*x = b
         * @param A d
         * @param x A template for cloning vectors needed within the solver.
         * @param b d
         * @param tol Tolerance (weighted 2-norm), iterative solvers only
         * @return error flag
         */
        int Solve(SUNMatrix A, N_Vector x, N_Vector b, realtype tol) const;

        /**
         * @brief Returns the last error flag encountered within the linear solver
         * @return error flag
         */
        [[nodiscard]] long int getLastFlag() const;

        /**
         * @brief Returns the integer and real workspace sizes for the linear solver
         * @param lenrwLS output argument for size of real workspace
         * @param leniwLS output argument for size of integer workspace
         * @return workspace size
         */
        int space(long int *lenrwLS, long int *leniwLS) const;

        /**
         * @brief Get the matrix A (matrix solvers only).
         * @return A
         */
        [[nodiscard]] virtual SUNMatrix getMatrix() const;

    protected:
        /**
         * @brief Performs linear solver initialization (assumes that all
         * solver-specific options have been set).
         * @return error code
         */
        int initialize();

        /** Wrapped solver */
        SUNLinearSolver solver_ {nullptr};
    };

    /**
     * @brief SUNDIALS band direct solver.
     */
    class SUNLinSolBand : public SUNLinSolWrapper {
    public:
        /**
         * @brief Create solver using existing matrix A without taking ownership of
         * A.
         * @param x A template for cloning vectors needed within the solver.
         * @param A square matrix
         */
        SUNLinSolBand(N_Vector x, SUNMatrix A);

        /**
         * @brief Create new band solver and matrix A.
         * @param x A template for cloning vectors needed within the solver.
         * @param ubw upper bandwidth of band matrix A
         * @param lbw lower bandwidth of band matrix A
         */
        SUNLinSolBand(Vector const &x, int ubw, int lbw);

        [[nodiscard]] SUNMatrix getMatrix() const override;

    private:
        /** Matrix A for solver, only if created by here. */
        SUNMatrixWrapper A_;
    };


    /**
     * @brief SUNDIALS dense direct solver.
     */
    class SUNLinSolDense : public SUNLinSolWrapper {
    public:
        /**
         * @brief Create dense solver
         * @param x A template for cloning vectors needed within the solver.
         */
        explicit SUNLinSolDense(Vector const &x);

        [[nodiscard]] SUNMatrix getMatrix() const override;

    private:
        /** Matrix A for solver, only if created by here. */
        SUNMatrixWrapper A_;
    };

    /**
     * @brief A RAII wrapper for SUNNonLinearSolver structs which solve the
     * nonlinear system F (y) = 0 or G(y) = y.
     */
    class SUNNonLinSolWrapper {
    public:
        /**
         * @brief SUNNonLinSolWrapper from existing SUNNonlinearSolver
         * @param sol d
         */
        explicit SUNNonLinSolWrapper(SUNNonlinearSolver sol);

        virtual ~SUNNonLinSolWrapper();

        /**
         * @brief Copy constructor
         * @param other
         */
        SUNNonLinSolWrapper(const SUNNonLinSolWrapper &other) = delete;

        /**
         * @brief Move constructor
         * @param other d
         */
        SUNNonLinSolWrapper(SUNNonLinSolWrapper &&other) noexcept;

        /**
         * @brief Copy assignment
         * @param other d
         * @return d
         */
        SUNNonLinSolWrapper &operator=(const SUNNonLinSolWrapper &other) = delete;

        /**
         * @brief Move assignment
         * @param other d
         * @return d
         */
        SUNNonLinSolWrapper &operator=(SUNNonLinSolWrapper &&other) noexcept;

        /**
         * @brief Get the wrapped SUNNonlinearSolver
         * @return SUNNonlinearSolver
         */
        [[nodiscard]] SUNNonlinearSolver get() const;

        /**
         * @brief Get type ID of the solver
         * @return d
         */
        [[nodiscard]] SUNNonlinearSolver_Type getType() const;

        /**
         * @brief Setup solver
         * @param y  the initial iteration passed to the nonlinear solver.
         * @param mem the sundials integrator memory structure.
         * @return d
         */
        int setup(N_Vector y, void *mem);

        /**
         * @brief Solve the nonlinear system F (y) = 0 or G(y) = y.
         * @param y0 the initial iterate for the nonlinear solve. This must remain
         * unchanged throughout the solution process.
         * @param y the solution to the nonlinear system
         * @param w the solution error weight vector used for computing weighted
         * error norms.
         * @param tol the requested solution tolerance in the weighted root-mean-
         * squared norm.
         * @param callLSetup a flag indicating that the integrator recommends for
         * the linear solver setup function to be called.
         * @param mem the sundials integrator memory structure.
         * @return d
         */
        int Solve(N_Vector y0, N_Vector y, N_Vector w, realtype tol,
                  bool callLSetup, void *mem);

        /**
         * @brief Set function to evaluate the nonlinear residual function F(y) = 0
         * or the fixed point function G(y) = y
         * @param SysFn d
         * @return d
         */
        int setSysFn(SUNNonlinSolSysFn SysFn);

        /**
         * @brief Set linear solver setup function.
         * @param SetupFn d
         * @return d
         */
        int setLSetupFn(SUNNonlinSolLSetupFn SetupFn);

        /**
         * @brief Set linear solver solve function.
         * @param SolveFn d
         * @return d
         */
        int setLSolveFn(SUNNonlinSolLSolveFn SolveFn);

        /**
         * @brief Set function to test for convergence
         * @param CTestFn d
         * @param ctest_data d
         * @return d
         */
        int setConvTestFn(SUNNonlinSolConvTestFn CTestFn, void* ctest_data);

        /**
         * @brief Set maximum number of non-linear iterations
         * @param maxiters d
         * @return d
         */
        int setMaxIters(int maxiters);

        /**
         * @brief getNumIters
         * @return d
         */
        [[nodiscard]] long int getNumIters() const;

        /**
         * @brief getCurIter
         * @return d
         */
        [[nodiscard]] int getCurIter() const;

        /**
         * @brief getNumConvFails
         * @return d
         */
        [[nodiscard]] long int getNumConvFails() const;

    protected:
        /**
         * @brief initialize
         */
        void initialize();

        /** the wrapper solver */
        SUNNonlinearSolver solver = nullptr;
    };


    /**
     * @brief SUNDIALS Newton non-linear solver to solve F (y) = 0.
     */
    class SUNNonLinSolNewton : public SUNNonLinSolWrapper {
    public:
        /**
         * @brief Create Newton solver
         * @param x A template for cloning vectors needed within the solver.
         */
        explicit SUNNonLinSolNewton(N_Vector x);

        /**
         * @brief Create Newton solver for enabled sensitivity analysis
         * @param count Number of vectors in the nonlinear solve. When integrating
         * a system containing Ns sensitivities the value of count is:
         *    - Ns+1 if using a simultaneous corrector approach.
         *    - Ns if using a staggered corrector approach.
         * @param x A template for cloning vectors needed within the solver.
         */
        SUNNonLinSolNewton(int count, N_Vector x);

        /**
         * @brief Get function to evaluate the nonlinear residual function F(y) = 0
         * @param SysFn d
         * @return d
         */
        int getSysFn(SUNNonlinSolSysFn *SysFn) const;
    };


    /**
     * @brief SUNDIALS Fixed point non-linear solver to solve G(y) = y.
     */
    class SUNNonLinSolFixedPoint : public SUNNonLinSolWrapper {
    public:
        /**
         * @brief Create fixed-point solver
         * @param x template for cloning vectors needed within the solver.
         * @param m number of acceleration vectors to use
         */
        explicit SUNNonLinSolFixedPoint(const_N_Vector x, int m = 0);

        /**
         * @brief Create fixed-point solver for use with sensitivity analysis
         * @param count Number of vectors in the nonlinear solve. When integrating
         * a system containing Ns sensitivities the value of count is:
         *    - Ns+1 if using a simultaneous corrector approach.
         *    - Ns if using a staggered corrector approach.
         * @param x template for cloning vectors needed within the solver.
         * @param m number of acceleration vectors to use
         */
        SUNNonLinSolFixedPoint(int count, const_N_Vector x, int m = 0);

        /**
         * @brief Get function to evaluate the fixed point function G(y) = y
         * @param SysFn s
         * @return s
         */
        int getSysFn(SUNNonlinSolSysFn *SysFn) const;
    };

} // namespace suneigen

#endif //SUNEIGEN_SUNDIALS_LINSOL_WRAPPER_H
