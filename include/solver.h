
#ifndef SUNEIGEN_SOLVER_H
#define SUNEIGEN_SOLVER_H

#include <utility>
#include <functional>
#include <memory>
#include <vector>
#include <chrono>
#include <sundials/sundials_types.h>

#include "defines.h"
#include "vector.h"
#include "sundials_linsol_wrapper.h"
#include "model.h"

namespace suneigen {

    class ReturnData;
    class Model;
    class Solver;

    class Solver {
    public:

        /** Type of what is passed to Sundials solvers as user_data */
        using user_data_type = std::pair<Model*, Solver const*>;

        Solver() = default;

        /**
         * @brief Constructor
         * @param app SunEigen application context
         */
        // explicit Solver(SunApplication* app);

        /**
         * @brief Default destructor
         */
        virtual ~Solver();

        /**
         * @brief Copy constructor
         * @param other other.
         */
        Solver(const Solver& other);

        /**
         * @brief Copy assignment
         * @param other other.
         */
        Solver& operator=(const Solver& other) = delete;

        /**
         * @brief Clone this instance
         * @return The clone
         */
        virtual Solver *clone() const = 0;

        /**
         * @brief runs a forward simulation until the specified timepoint
         *
         * @param tout next timepoint
         * @return status flag
         */
        int run(realtype tout) const;

        /**
         * @brief Initializes the ami memory object and applies specified options
         * @param t0 initial timepoint
         * @param model pointer to the model instance
         * @param x0 initial states
         * @param dx0 initial derivative states
         * @param sx0 initial state sensitivities
         * @param sdx0 initial derivative state sensitivities
         */
         void setup(realtype t0, Model *model, const Vector& x0,
                 const Vector& dx0, const VectorArray& sx0,
                 const VectorArray& sdx0) const;

        /**
        * getRootInfo extracts information which event occurred
        *
        * @param rootsfound array with flags indicating whether the respective
        * event occurred
        */
        virtual void getRootInfo(int* rootsfound) const = 0;

        /**
        * @brief Calculates consistent initial conditions, assumes initial
        * states to be correct (DAE only)
        *
        * @param tout1 next timepoint to be computed (sets timescale)
        */
        virtual void calcIC(realtype tout1) const = 0;

        /**
         * @brief Disable rootfinding
         */
        virtual void turnOffRootFinding() const = 0;


        /**
        * @brief Return current sensitivity method
        * @return method enum
        */
        SensitivityMethod getSensitivityMethod() const;

        /**
         * @brief Set sensitivity method
         * @param sensi_meth g
         */
        void setSensitivityMethod(SensitivityMethod sensi_meth);

        /**
         * @brief Get maximum number of allowed linear steps per Newton step for
         * steady state computation
         * @return dsf
         */
        size_t getNewtonMaxLinearSteps() const;

        /**
         * @brief Set maximum number of allowed linear steps per Newton step for
         * steady state computation
         * @param newton_maxlinsteps fdsfs
         */
        void setNewtonMaxLinearSteps(size_t newton_maxlinsteps);

        /**
         * @brief Get a state of the damping factor used in the Newton solver
         * @return ds
         */
        NewtonDampingFactorMode getNewtonDampingFactorMode() const;

        /**
         * @brief Turn on/off a damping factor in the Newton method
         * @param dampingFactorMode sdfs
         */
        void setNewtonDampingFactorMode(NewtonDampingFactorMode dampingFactorMode);

        /**
         * @brief Get a lower bound of the damping factor used in the Newton solver
         * @return dsf
         */
        realtype getNewtonDampingFactorLowerBound() const;

        /**
         * @brief Set a lower bound of the damping factor in the Newton solver
         * @param dampingFactorLowerBound sdfs
         */
        void setNewtonDampingFactorLowerBound(realtype dampingFactorLowerBound);

        /**
         * @brief Get sensitivity order
         * @return sensitivity order
         */
        SensitivityOrder getSensitivityOrder() const;

        /**
         * @brief Set the sensitivity order
         * @param sensi sensitivity order
         */
        void setSensitivityOrder(SensitivityOrder sensi);

        /**
         * @brief Get the relative tolerances for the forward problem
         *
         * Same tolerance is used for the backward problem if not specified
         * differently via setRelativeToleranceASA.
         *
         * @return relative tolerances
         */
        realtype getRelativeTolerance() const;

        /**
         * @brief Sets the relative tolerances for the forward problem
         *
         * Same tolerance is used for the backward problem if not specified
         * differently via setRelativeToleranceASA.
         *
         * @param rtol relative tolerance (non-negative number)
         */
        void setRelativeTolerance(realtype rtol);

        /**
         * @brief Get the absolute tolerances for the forward problem
         *
         * Same tolerance is used for the backward problem if not specified
         * differently via setAbsoluteToleranceASA.
         *
         * @return absolute tolerances
         */
        realtype getAbsoluteTolerance() const;

        /**
         * @brief Sets the absolute tolerances for the forward problem
         *
         * Same tolerance is used for the backward problem if not specified
         * differently via setAbsoluteToleranceASA.
         *
         * @param atol absolute tolerance (non-negative number)
         */
        void setAbsoluteTolerance(realtype atol);

        /**
         * @brief Returns the relative tolerances for the forward sensitivity analysis
         * problem
         * @return relative tolerances
         */
        realtype getRelativeToleranceFSA() const;

        /**
         * @brief Sets the relative tolerances for the forward sensitivity analysis
         * @param rtol relative tolerance (non-negative number)
         */
        void setRelativeToleranceFSA(realtype rtol);

        /**
         * @brief Returns the absolute tolerances for the forward sensitivity analysis
         * @return absolute tolerances
         */
        realtype getAbsoluteToleranceFSA() const;

        /**
         * @brief Sets the absolute tolerances for the forward sensitivity analysis
         * @param atol absolute tolerance (non-negative number)
         */
        void setAbsoluteToleranceFSA(realtype atol);

        /**
         * @brief returns the maximum number of solver steps for the forward
         * problem
         * @return maximum number of solver steps
         */
        size_t getMaxSteps() const;

        /**
         * @brief sets the maximum number of solver steps for the forward problem
         * @param maxsteps maximum number of solver steps (positive number)
         */
        void setMaxSteps(size_t maxsteps);

        /**
         * @brief Returns the maximum time allowed for integration
         * @return Time in seconds
         */
        double getMaxTime() const;

        /**
         * @brief Set the maximum time allowed for integration
         * @param maxtime Time in seconds
         */
        void setMaxTime(double maxtime);

        /**
         * @brief Start timer for tracking integration time
        */
        void startTimer() const;

        /**
         * @brief Check whether maximum integration time was exceeded
         * @return True if the maximum integration time was exceeded,
         * false otherwise.
         */
        bool timeExceeded() const;

        /**
         * @brief returns the linear system multistep method
         * @return linear system multistep method
         */
        LinearMultistepMethod getLinearMultistepMethod() const;

        /**
         * @brief sets the linear system multistep method
         * @param lmm linear system multistep method
         */
        void setLinearMultistepMethod(LinearMultistepMethod lmm);

        /**
         * @brief returns the nonlinear system solution method
         * @return iteration counter
         */
        NonlinearSolverIteration getNonlinearSolverIteration() const;

        /**
         * @brief sets the nonlinear system solution method
         * @param iter nonlinear system solution method
         */
        void setNonlinearSolverIteration(NonlinearSolverIteration iter);

        /**
         * @brief getInterpolationType
         * @return interp type
         */
        InterpolationType getInterpolationType() const;

        /**
         * @brief sets the interpolation of the forward solution that is used for
         * the backwards problem
         * @param interpType interpolation type
         */
        void setInterpolationType(InterpolationType interpType);

        /**
         * @brief getLinearSolver
         * @return linear solver
         */
        LinearSolver getLinearSolver() const;

        /**
         * @brief setLinearSolver
         * @param linsol ff
         */
        void setLinearSolver(LinearSolver linsol);

        /**
         * @brief returns the internal sensitivity method
         * @return internal sensitivity method
         */
        InternalSensitivityMethod getInternalSensitivityMethod() const;

        /**
         * @brief sets the internal sensitivity method
         * @param ism internal sensitivity method
         */
        void setInternalSensitivityMethod(InternalSensitivityMethod ism);

        /**
         * @brief write solution from forward simulation
         * @param t time
         * @param x state
         * @param dx derivative state
         * @param sx state sensitivity
         * @param xQ quadrature
         */
        void writeSolution(realtype *t, Vector &x, Vector &dx,
                           VectorArray &sx, Vector &xQ) const;

        /**
         * @brief Access state solution at time t
         * @param t time
         * @return x or interpolated solution dky
         */
        const Vector& getState(realtype t) const;

        /**
         * @brief Access derivative state solution at time t
         * @param t time
         * @return dx or interpolated solution dky
         */
        const Vector &getDerivativeState(realtype t) const;

        /**
         * @brief Access state sensitivity solution at time t
         * @param t time
         * @return (interpolated) solution sx
         */
        const VectorArray &getStateSensitivity(realtype t) const;

        /**
         * @brief Access quadrature solution at time t
         * @param t time
         * @return (interpolated) solution xQ
         */
        const Vector& getQuadrature(realtype t) const;

        /**
         * @brief Reinitializes the states in the solver after an event occurrence
         *
         * @param t0 reinitialization timepoint
         * @param yy0 initial state variables
         * @param yp0 initial derivative state variables (DAE only)
         */
        virtual void reInit(realtype t0, const Vector &yy0,
                            const Vector &yp0) const = 0;

        /**
         * @brief Reinitializes the state sensitivities in the solver after an
         * event occurrence
         *
         * @param yyS0 new state sensitivity
         * @param ypS0 new derivative state sensitivities (DAE only)
         */
        virtual void sensReInit(const VectorArray &yyS0,
                                const VectorArray &ypS0) const = 0;

        /**
         * @brief Switches off computation of  state sensitivities without
         * deallocating the memory for sensitivities
         */
        virtual void sensToggleOff() const = 0;

        /**
         * @brief current solver timepoint
         * @return t
         */
        realtype gett() const;

        /**
         * @brief Reads out the CPU time needed for forward solve
         * @return cpu_time
         */
        realtype getCpuTime() const;

        /**
         * @brief number of states with which the solver was initialized
         * @return x.getLength()
         */
        size_t nx() const;

        /**
         * @brief number of parameters with which the solver was initialized
         * @return sx.getLength()
         */
        unsigned int np() const;

        /**
         * @brief number of quadratures with which the solver was initialized
         * @return xQB.getLength()
         */
        size_t nquad() const;

        /**
         * @brief check if FSA is being computed
         * @return flag
         */
        bool computingFSA() const {
            return getSensitivityOrder() >= SensitivityOrder::first &&
                   getSensitivityMethod() == SensitivityMethod::forward && np() > 0;
        }

        /**
         * @brief check if ASA is being computed
         * @return flag
         */
        bool computingASA() const {
            return getSensitivityOrder() >= SensitivityOrder::first &&
                   getSensitivityMethod() == SensitivityMethod::adjoint && np() > 0;
        }

        /**
         * @brief Resets vectors containing diagnosis information
         */
        void resetDiagnosis() const;

        /**
         * @brief Stores diagnosis information from solver memory block for forward problem
         */
        void storeDiagnosis() const;

        /**
         * @brief Accessor ns
         * @return ns
         */
        std::vector<size_t> const& getNumSteps() const {
            return ns_;
        }

        /**
         * @brief Accessor nrhs
         * @return nrhs
         */
        std::vector<size_t> const& getNumRhsEvals() const {
            return nrhs_;
        }

        /**
         * @brief Accessor netf
         * @return netf
         */
        std::vector<size_t> const& getNumErrTestFails() const {
            return netf_;
        }

        /**
         * @brief Accessor netfB
         * @return netfB
         */
        std::vector<size_t> const& getNumErrTestFailsB() const {
            return netfB_;
        }

        /**
         * @brief Accessor nnlscf
         * @return nnlscf
         */
        std::vector<size_t> const& getNumNonlinSolvConvFails() const {
            return nnlscf_;
        }

        /**
         * @brief Accessor nnlsi
         * @return nnlscf
         */
        std::vector<size_t> const& getNumNonlinSolvIters() const {
            return nnlsi_;
        }

        /**
         * @brief Accessor nje
         * @return nje
         */
        std::vector<size_t> const& getNumJacEvals() const {
            return nje_;
        }

        /**
         * @brief Accessor nnlscfB
         * @return nnlscfB
         */
        std::vector<size_t> const& getNumNonlinSolvConvFailsB() const {
            return nnlscfB_;
        }

        /**
         * @brief Accessor order
         * @return order
         */
        std::vector<size_t> const& getLastOrder() const {
            return order_;
        }

        /** maximum number of allowed linear steps per Newton step for steady state
        * computation */
        size_t newton_maxlinsteps_ {0L};

        /**
         * @brief Get the solver type used (cvodes or idas)
         * @return A string representing the solver used.
         */
        virtual std::string getSolverType() const = 0;

        /**
         * @brief Check equality of data members excluding solver memory
         * @param a Solver
         * @param b Solver
         * @return success or not
         */
        friend bool operator==(const Solver &a, const Solver &b);

    protected:

        /**
         * @brief Sets a timepoint at which the simulation will be stopped
         *
         * @param tstop timepoint until which simulation should be performed
         */
        virtual void setStopTime(realtype tstop) const = 0;

        /**
         * @brief Solves the forward problem until a predefined timepoint
         *
         * @param tout timepoint until which simulation should be performed
         * @param itask task identifier, can be CV_NORMAL or CV_ONE_STEP
         * @return status flag indicating success of execution
         */
        virtual int solve(realtype tout, int itask) const = 0;

        /**
         * @brief Solves the forward problem until a predefined timepoint
         * (adjoint only)
         *
         * @param tout timepoint until which simulation should be performed
         * @param itask task identifier, can be CV_NORMAL or CV_ONE_STEP
         * @param ncheckPtr pointer to a number that counts the internal
         * checkpoints
         * @return status flag indicating success of execution
         */
        // virtual int solveF(realtype tout, int itask, int *ncheckPtr) const = 0;

        /**
         * @brief reInitPostProcessF postprocessing of the solver memory after a
         * discontinuity in the forward problem
         * @param tnext next timepoint (defines integration direction)
         */
        virtual void reInitPostProcessF(realtype tnext) const = 0;

        /**
         * @brief extracts the state sensitivity at the current timepoint from
         * solver memory and writes it to the sx member variable
         */
        virtual void getSens() const = 0;

        /**
         * @brief extracts the quadrature at the current timepoint from solver
         * memory and writes it to the xQ member variable
         *
         * @param t timepoint for quadrature extraction
         */
        virtual void getQuad(realtype &t) const = 0;

        /**
         * @brief Initializes the states at the specified initial timepoint
         *
         * @param t0 initial timepoint
         * @param x0 initial states
         * @param dx0 initial derivative states
         */
        virtual void init(realtype t0, const Vector &x0,
                          const Vector &dx0) const = 0;

        /**
         * @brief Initializes the rootfinding for events
         *
         * @param ne number of different events
         */
        virtual void rootInit(int ne) const = 0;

        /**
         * @brief Set the dense Jacobian function
         */
        virtual void setDenseJacFn() const = 0;

        /**
         * @brief sets the sparse Jacobian function
         */
        virtual void setSparseJacFn() const = 0;

        /**
         * @brief Create specifies solver method and initializes solver memory for
         * the forward problem
         */
        virtual void allocateSolver() const = 0;

        /**
         * @brief sets scalar relative and absolute tolerances for the forward
         * problem
         *
         * @param rtol relative tolerances
         * @param atol absolute tolerances
         */
        virtual void setSStolerances(double rtol, double atol) const = 0;

        /**
         * @brief activates sets scalar relative and absolute tolerances for the
         * sensitivity variables
         *
         * @param rtol relative tolerances
         * @param atol array of absolute tolerances for every sensitivity variable
         */
        virtual void setSensSStolerances(double rtol, const double *atol) const = 0;

        /**
         * SetSensErrCon specifies whether error control is also enforced for
         * sensitivities for the forward problem
         *
         * @param error_corr activation flag
         */
        virtual void setSensErrCon(bool error_corr) const = 0;

        /**
         * @brief Attaches the error handler function (errMsgIdAndTxt)
         * to the solver
         *
         */
        virtual void setErrHandlerFn() const = 0;

        /**
         * @brief Attaches the user data to the forward problem
        */
        virtual void setUserData() const = 0;

        /**
        * @brief specifies the maximum number of steps for the forward problem
         *
         * @param mxsteps number of steps
         * @note in contrast to the SUNDIALS method, this sets the overall maximum, not the maximum between output times.
         */
        virtual void setMaxNumSteps(size_t mxsteps) const = 0;

        /**
         * @brief activates stability limit detection for the forward
         * problem
         *
         * @param stldet flag for stability limit detection (TRUE or FALSE)
         *
         */
        virtual void setStabLimDet(int stldet) const = 0;

        /**
         * @brief specify algebraic/differential components (DAE only)
         *
         * @param model model specification
         */
        virtual void setId(const Model *model) const = 0;

        /**
         * @brief deactivates error control for algebraic components (DAE only)
         *
         * @param flag deactivation flag
         */
        virtual void setSuppressAlg(bool flag) const = 0;

        /**
         * @brief specifies the scaling and indexes for sensitivity
         * computation
         *
         * @param p parameters
         * @param pbar parameter scaling constants
         * @param plist parameter index list
         */
        virtual void setSensParams(const realtype *p, const realtype *pbar,
                                   const int *plist) const = 0;

        /**
         * @brief interpolates the (derivative of the) solution at the requested
         * timepoint
         *
         * @param t timepoint
         * @param k derivative order
         */
        virtual void getDky(realtype t, int k) const = 0;

        /**
         * @brief interpolates the (derivative of the) solution at the requested
         * timepoint
         *
         * @param t timepoint
         * @param k derivative order
         */
        virtual void getSensDky(realtype t, int k) const = 0;

        /**
         * @brief interpolates the (derivative of the) solution at the requested
         * timepoint
         *
         * @param t timepoint
         * @param k derivative order
         */
        virtual void getQuadDky(realtype t, int k) const = 0;

        /**
         * @brief reports the number of solver steps
         *
         * @param sun_mem pointer to the solver memory instance (can be from
         * forward or backward problem)
         * @param numsteps output array
         */
        virtual void getNumSteps(const void *sun_mem, size_t *numsteps) const = 0;

        /**
         * @brief reports the number of nonlinear solver iterations.
         * @param sun_mem pointer to the solver memory instance (can be from
         * forward or backward problem)
         * @param numnonlinsolviters output array
         */
        virtual void getNumNonlinSolvIters(const void* sun_mem, size_t *numnonlinsolviters) const = 0;

        /**
         * @brief reports the number of jacobian evaluations.
         * @param sun_mem pointer to the solver memory instance (can be from
         * forward or backward problem)
         * @param numjacevals output array
         */
        virtual void getNumJacEvals(const void* sun_mem, size_t *numjacevals) const = 0;

        /**
         * @brief reports the number of right hand evaluations
         *
         * @param sun_mem pointer to the solver memory instance (can be from
         * forward or backward problem)
         * @param numrhsevals output array
         */
        virtual void getNumRhsEvals(const void *sun_mem,
                                    size_t *numrhsevals) const = 0;

        /**
         * @brief reports the number of local error test failures
         *
         * @param sun_mem pointer to the solver memory instance (can be from
         * forward or backward problem)
         * @param numerrtestfails output array
         */
        virtual void getNumErrTestFails(const void *sun_mem,
                                        size_t *numerrtestfails) const = 0;

        /**
         * @brief reports the number of nonlinear convergence failures
         *
         * @param sun_mem pointer to the solver memory instance (can be from
         * forward or backward problem)
         * @param numnonlinsolvconvfails output array
         */
        virtual void
        getNumNonlinSolvConvFails(const void *sun_mem,
                                  size_t *numnonlinsolvconvfails) const = 0;

        /**
         * @brief Reports the order of the integration method during the
         * last internal step
         *
         * @param sun_mem pointer to the solver memory instance (can be from
         * forward or backward problem)
         * @param order output array
         */
        virtual void getLastOrder(const void *sun_mem, size_t *order) const = 0;

        /**
         * @brief Initializes and sets the linear solver for the forward problem
         *
         * @param model pointer to the model object
         */
        void initializeLinearSolver(const Model *model) const;

        /**
         * @brief Sets the non-linear solver
         */
        void initializeNonLinearSolver() const;

        /**
         * @brief Sets the linear solver for the forward problem
         */
        virtual void setLinearSolver() const = 0;

        /**
         * @brief Set the non-linear solver for the forward problem
         */
        virtual void setNonLinearSolver() const = 0;

        /**
         * @brief Set the non-linear solver for sensitivities
         */
        virtual void setNonLinearSolverSens() const = 0;

        /**
         * Accessor function to the model stored in the user data
         *
         * @return user data model
         */
        virtual const Model *getModel() const = 0;

        /**
         * @brief checks whether memory for the forward problem has been allocated
         *
         * @return proxy for solverMemory->(cv|ida)_MallocDone
         */
        bool getInitDone() const;

        /**
         * @brief checks whether memory for forward sensitivities has been allocated
         *
         * @return proxy for solverMemory->(cv|ida)_SensMallocDone
         */
        bool getSensInitDone() const;

        /**
         * @brief checks whether memory for forward interpolation has been allocated
         *
         * @return proxy for solverMemory->(cv|ida)_adjMallocDone
         */
        bool getAdjInitDone() const;

        /**
         * @brief resets solverMemory and solverMemoryB
         * @param nx new number of state variables
         * @param nplist new number of sensitivity parameters
         * @param nquad new number of quadratures (only differs from nplist for
         * higher order sensitivity computation)
         */
        void resetMutableMemory(size_t nx, size_t nplist, size_t nquad) const;

        /**
         * @brief updates solver tolerances according to the currently specified
         * member variables
         */
        void applyTolerances() const;

        /**
         * @brief updates FSA solver tolerances according to the currently
         * specified member variables
         */
        void applyTolerancesFSA() const;

        /**
         * @brief updates all sensitivity solver tolerances according to the
         * currently specified member variables
         */
        void applySensitivityTolerances() const;

        /** pointer to solver memory block */
        mutable std::unique_ptr<void, std::function<void(void *)>> solver_memory_;

        /** Sundials user_data */
        mutable user_data_type user_data;

        /** internal sensitivity method flag used to select the sensitivity solution
         * method. Only applies for Forward Sensitivities. */
        InternalSensitivityMethod ism_ {InternalSensitivityMethod::simultaneous};

        /** specifies the linear multistep method.
         */
        LinearMultistepMethod lmm_ {LinearMultistepMethod::BDF};

        /**
         * specifies the type of nonlinear solver iteration
         */
        NonlinearSolverIteration iter_ {NonlinearSolverIteration::newton};

        /** interpolation type for the forward problem solution which
         * is then used for the backwards problem.
        */
        InterpolationType interp_type_ {InterpolationType::hermite};

        /** maximum number of allowed integration steps */
        size_t maxsteps_ {10000};

        /** Maximum wall-time for integration in seconds */
        std::chrono::duration<double, std::ratio<1>> maxtime_ {std::chrono::duration<double>::max()};

        /** Time at which solver timer was started */
        mutable std::chrono::time_point<std::chrono::system_clock> starttime_;

        /** linear solver for the forward problem */
        mutable std::unique_ptr<SUNLinSolWrapper> linear_solver_;

        /** non-linear solver for the forward problem */
        mutable std::unique_ptr<SUNNonLinSolWrapper> non_linear_solver_;

        /** non-linear solver for the sensitivities */
        mutable std::unique_ptr<SUNNonLinSolWrapper> non_linear_solver_sens_;

        /** flag indicating whether the forward solver has been called */
        mutable bool solver_was_called_F_ {false};

        /** flag indicating whether the backward solver has been called */
        mutable bool solver_was_called_B_ {false};

        /**
         * @brief sets that memory for the forward problem has been allocated
         */
        void setInitDone() const;

        /**
         * @brief sets that memory for forward sensitivities has been allocated
         */
        void setSensInitDone() const;

        /**
         * @brief Sets sensitivity method
         * @param sensi_meth new value for sensi_meth
         */
        void checkSensitivityMethod(SensitivityMethod sensi_meth) const;

        /** integration time of the forward problem */
        mutable realtype t_ {std::nan("")};

        /** state (dimension: nx_solver) */
        mutable Vector x_ {0};

        /** state interface variable (dimension: nx_solver) */
        mutable Vector dky_ {0};

        /** state derivative dummy (dimension: nx_solver) */
        mutable Vector dx_ {0};

        /** state sensitivities interface variable (dimension: nx_solver x nplist) */
        mutable VectorArray sx_ {0, 0};

        /** state derivative sensitivities dummy (dimension: nx_solver x nplist) */
        mutable VectorArray sdx_ {0, 0};

        /** forward quadrature interface variable (dimension: nx_solver) */
        mutable Vector xQ_ {0};

        /** adjoint quadrature interface variable (dimension: nJ x nplist) */
        mutable Vector xQB_ {0};

        /** flag to force reInitPostProcessF before next call to solve */
        mutable bool force_reinit_postprocess_F_ {false};

    private:

        /**
         * @brief applies total number of steps for next solver call
         */
        void apply_max_num_steps() const;

        /** method for sensitivity computation */
        SensitivityMethod sensi_meth_ {SensitivityMethod::forward};

        /** flag controlling stability limit detection */
        booleantype stldet_ {true};

        /** linear solver specification */
        LinearSolver linsol_ {LinearSolver::SuperLU};

        /** maximum number of allowed Newton steps for steady state computation */
        // long int newton_maxsteps_ {0L};

        /** Damping factor state used int the Newton method */
        NewtonDampingFactorMode newton_damping_factor_mode_{NewtonDampingFactorMode::on};

        /** CPU time, forward solve */
        mutable realtype cpu_time_ {0.0};

        /** absolute tolerances for forward sensitivity integration */
        realtype atol_fsa_ {std::nan("")};

        /** relative tolerances for forward sensitivity integration */
        realtype rtol_fsa_ {std::nan("")};

        /** Lower bound of the damping factor. */
        realtype newton_damping_factor_lower_bound_ {1e-8};

        /** absolute tolerances for integration */
        realtype atol_ {1e-16};

        /** relative tolerances for integration */
        realtype rtol_ {1e-8};

        /** number of integration steps forward problem (dimension: nt) */
        mutable std::vector<size_t> ns_;

        /** number of integration steps backward problem (dimension: nt) */
        mutable std::vector<size_t> nsB_;

        /** number of right hand side evaluations forward problem (dimension: nt) */
        mutable std::vector<size_t> nrhs_;

        /** number of right hand side evaluations backward problem (dimension: nt) */
        mutable std::vector<size_t> nrhsB_;

        /** number of error test failures forward problem (dimension: nt) */
        mutable std::vector<size_t> netf_;

        /** number of error test failures backward problem (dimension: nt) */
        mutable std::vector<size_t> netfB_;

        /**
         * number of non-linear solver convergence failures forward problem (dimension:
         * nt) */
        mutable std::vector<size_t> nnlscf_;

        /**
         * Number of non-linear solver iterations
         */
        mutable std::vector<size_t> nnlsi_;

        /**
         * Number of jacobian evaluations.
         */
        mutable std::vector<size_t> nje_;
        /**
         * number of non-linear solver convergence failures backward problem (dimension:
         * nt) */
        mutable std::vector<size_t> nnlscfB_;

        /** employed order forward problem (dimension: nt) */
        mutable std::vector<size_t> order_;

        /** flag indicating whether sensitivities are supposed to be computed */
        SensitivityOrder sensi_ {SensitivityOrder::none};

        /** flag indicating whether init was called */
        mutable bool initialized_ {false};

        /** flag indicating whether sensInit1 was called */
        mutable bool sens_initialized_ {false};

        /** flag indicating whether adjInit was called */
        mutable bool adj_initialized_ {false};

        /** flag indicating whether (forward) quadInit was called */
        mutable bool quad_initialized_ {false};

    };

    /**
     * @brief Extracts diagnosis information from solver memory block and
     * passes them to the specified output function
     *
     * @param error_code error identifier
     * @param module name of the module in which the error occurred
     * @param function name of the function in which the error occurred
     * @param msg error message
     * @param eh_data amici::Solver as void*
     */
    void wrapErrHandlerFn(int error_code, const char *module, const char *function,
                          char *msg, void *eh_data);

} // namespace suneigen

#endif //SUNEIGEN_SOLVER_H
