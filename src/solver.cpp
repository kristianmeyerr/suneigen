//
// Created by Kristian Meyer on 04/11/2021.
//

#include <sundials_linsol_wrapper.h>
#include "solver.h"
#include "exception.h"
#include "suneigen.h"

#include <ctime>

namespace suneigen {

    //Solver::Solver(SunApplication* app_) : app(app_){}

    void Solver::apply_max_num_steps() const {
        // set remaining steps, setMaxNumSteps only applies to a single call of solve
        size_t cursteps;
        getNumSteps(solver_memory_.get(), &cursteps);
        if (maxsteps_ <= cursteps)
            throw SunException("Reached maximum number of steps %ld before reaching "
                               "tout at t=%g.", maxsteps_, t_);
        setMaxNumSteps(maxsteps_ - cursteps);
    }

    void Solver::apply_max_error_test_fails() const {
        setMaxErrTestFails(maxErrorTestFails_);
    }

    Solver::~Solver() = default;

    Solver::Solver(const Solver& other) : atol_(other.atol_) {}

    int Solver::run(const realtype tout) const {
        setStopTime(tout);
        clock_t starttime = clock();

        int status = SUNEIGEN_SUCCESS;

        apply_max_num_steps();
        if (nx() > 0) {
            if (getAdjInitDone()) {
                // Solve
            } else {
                status = solve(tout, SUNEIGEN_NORMAL);
            }
        } else {
            t_ = tout;
        }
        cpu_time_ += static_cast<realtype>(((clock() - starttime) * 1000)) / CLOCKS_PER_SEC;
        return status;
    }

    realtype Solver::getRelativeTolerance() const {
        return rtol_;
    }

    void Solver::setRelativeTolerance(const realtype rtol) {
        if (rtol < 0)
            throw SunException("rtol must be a non-negative number");

        rtol_ = rtol;

        if (getInitDone()) {
            applyTolerances();
            applySensitivityTolerances();
        }
    }

    realtype Solver::getAbsoluteTolerance() const {
        return atol_;
    }

    void Solver::setAbsoluteTolerance(realtype atol) {
        if (atol < 0)
            throw SunException("atol must be a non-negative number");

        atol_ = atol;

        if (getInitDone()) {
            applyTolerances();
            applySensitivityTolerances();
        }
    }

    size_t Solver::getMaxSteps() const { return maxsteps_; }

    void Solver::setMaxSteps(const size_t maxsteps) {
        if (maxsteps <= 0)
            throw SunException("maxsteps must be a positive number");

        maxsteps_ = maxsteps;

        if (getAdjInitDone())
            resetMutableMemory(nx(), np(), nquad());

    }

    void Solver::setMaxErrorTests(size_t maxfails) {
        if (maxfails <= 0)
            throw SunException("maxsteps must be a positive number");

        maxErrorTestFails_ = maxfails;

        if (getAdjInitDone())
            resetMutableMemory(nx(), np(), nquad());

    }

    bool Solver::getAdjInitDone() const { return adj_initialized_; }

    double Solver::getMaxTime() const { return maxtime_.count(); }

    void Solver::setMaxTime(double maxtime)
    {
        maxtime_ = std::chrono::duration<double>(maxtime);
    }

    void Solver::startTimer() const
    {
        starttime_ = std::chrono::system_clock::now();
    }

    bool Solver::timeExceeded() const
    {
        return std::chrono::system_clock::now() - starttime_ > maxtime_;
    }

    LinearMultistepMethod Solver::getLinearMultistepMethod() const { return lmm_; }

    void Solver::setLinearMultistepMethod(const LinearMultistepMethod lmm) {
        if (solver_memory_)
            resetMutableMemory(nx(), np(), nquad());
        lmm_ = lmm;
    }

    NonlinearSolverIteration Solver::getNonlinearSolverIteration() const {
        return iter_;
    }

    void Solver::setNonlinearSolverIteration(const NonlinearSolverIteration iter) {
        if (solver_memory_)
            resetMutableMemory(nx(), np(), nquad());
        iter_ = iter;
    }

    InterpolationType Solver::getInterpolationType() const { return interp_type_; }

    void Solver::setInterpolationType(const InterpolationType interpType) {

        /*
        if (!solver_memory_B_.empty())
            resetMutableMemory(nx(), nplist(), nquad());
        */
        interp_type_ = interpType;
    }

    LinearSolver Solver::getLinearSolver() const { return linsol_; }

    void Solver::setLinearSolver(LinearSolver linsol) {
        if (solver_memory_)
            resetMutableMemory(nx(), np(), nquad());
        linsol_ = linsol;
    }

    InternalSensitivityMethod Solver::getInternalSensitivityMethod() const {
        return ism_;
    }

    size_t Solver::getNewtonMaxLinearSteps() const { return newton_maxlinsteps_; }

    void Solver::setNewtonMaxLinearSteps(const size_t newton_maxlinsteps) {
        newton_maxlinsteps_ = newton_maxlinsteps;
    }

    NewtonDampingFactorMode Solver::getNewtonDampingFactorMode() const { return newton_damping_factor_mode_; }

    void Solver::setNewtonDampingFactorMode(NewtonDampingFactorMode dampingFactorMode) {
        newton_damping_factor_mode_ = dampingFactorMode;
    }

    realtype Solver::getNewtonDampingFactorLowerBound() const { return newton_damping_factor_lower_bound_; }

    void Solver::setNewtonDampingFactorLowerBound(realtype dampingFactorLowerBound) {
        newton_damping_factor_lower_bound_ = dampingFactorLowerBound;
    }

    realtype Solver::getCpuTime() const {
        return cpu_time_;
    }

    void Solver::setInternalSensitivityMethod(const InternalSensitivityMethod ism) {
        if (solver_memory_)
            resetMutableMemory(nx(), np(), nquad());
        ism_ = ism;
    }


    void Solver::resetMutableMemory(const size_t nx, const size_t np,
                                    const size_t nquad) const {
        solver_memory_ = nullptr;
        initialized_ = false;
        adj_initialized_ = false;
        sens_initialized_ = false;
        quad_initialized_ = false;
        solver_was_called_F_ = false;
        solver_was_called_B_ = false;

        x_ = Vector(nx);  // @todo Is the copy or move assignment operator called here?
        dx_ = Vector(nx);
        sx_ = VectorArray(nx, np);
        sdx_ = VectorArray(nx, np);

        //xB_ = Vector(nx);
        //dxB_ = Vector(nx);
        xQB_ = Vector(nquad);
        xQ_ = Vector(nx);

        //solver_memory_B_.clear();
        //initializedB_.clear();
        //initializedQB_.clear();
    }

    void Solver::setup(realtype t0, Model *model, const Vector& x0,
                       const Vector& dx0, const VectorArray& sx0,
                       const VectorArray& sdx0) const {
        (void)sx0; (void)sdx0;

        if (nx() != model->nx || np() != model->np())
            resetMutableMemory(model->nx, model->np(), model->np());

        /* Create solver memory object if necessary */
        allocateSolver();
        if (!solver_memory_)
            throw SunException("Failed to allocated solver memory!");

        /* Initialize CVodes/IDAs solver*/
        init(t0, x0, dx0);

        /* Clear diagnosis storage */
        resetDiagnosis();

        /* Apply stored tolerances to sundials solver */
        applyTolerances();

        /* Set optional inputs */
        setErrHandlerFn();

        /* Attaches userdata */
        user_data = std::make_pair(model, this);
        setUserData();

        /* activates stability limit detection */
        setStabLimDet(stldet_);

        rootInit(static_cast<int>(model->ne));

        if (nx() == 0)
            return;

        initializeLinearSolver(model);
        initializeNonLinearSolver();

        // Work on sensitivities

        setId(model);
        setSuppressAlg(true);
        /* calculate consistent DAE initial conditions (no effect for ODE) */
        if (model->nt() > 1)
            calcIC(model->getTimepoint(1));
    }

    void Solver::writeSolution(realtype* t, Vector &x, Vector &dx,
                               VectorArray& sx, Vector &xQ) const {
        *t = gett();
        if (quad_initialized_)
            xQ.copy(getQuadrature(*t));
        if (sens_initialized_)
            sx.copy(getStateSensitivity(*t));
        x.copy(getState(*t));
        dx.copy(getDerivativeState(*t));
    }

    const Vector &Solver::getQuadrature(realtype t) const {
        if (quad_initialized_) {
            if (solver_was_called_F_) {
                if (t == t_) {
                    getQuad(t);
                    return xQ_;
                }
                getQuadDky(t, 0);
            }
        } else {
            xQ_.zero();
        }
        return xQ_;
    }

    realtype Solver::gett() const { return t_; }

    unsigned int Solver::np() const { return static_cast<unsigned int>(sx_.getLength()); }

    size_t Solver::nx() const { return x_.getLength(); }

    size_t Solver::nquad() const { return xQB_.getLength(); }

    bool Solver::getInitDone() const { return initialized_; }

    void Solver::setInitDone() const { initialized_ = true; }

    void Solver::resetDiagnosis() const {
        ns_.clear();
        nrhs_.clear();
        netf_.clear();
        nnlscf_.clear();
        nnlsi_.clear();
        nje_.clear();
        order_.clear();

        nsB_.clear();
        nrhsB_.clear();
        netfB_.clear();
        nnlscfB_.clear();
    }

    bool operator==(const Solver &a, const Solver &b) {
        bool is_equal = true;
        if (typeid(a) != typeid(b))
            is_equal = false;
        return is_equal;
    }

    void Solver::storeDiagnosis() const {
        if (!solver_was_called_F_ || !solver_memory_) {
            ns_.push_back(0);
            nrhs_.push_back(0);
            netf_.push_back(0);
            nnlscf_.push_back(0);
            nnlsi_.push_back(0);
            nje_.push_back(0);
            order_.push_back(0);
            return;
        }

        size_t lnumber;
        getNumSteps(solver_memory_.get(), &lnumber);
        ns_.push_back(lnumber);

        getNumRhsEvals(solver_memory_.get(), &lnumber);
        nrhs_.push_back(lnumber);

        getNumErrTestFails(solver_memory_.get(), &lnumber);
        netf_.push_back(lnumber);

        getNumNonlinSolvConvFails(solver_memory_.get(), &lnumber);
        nnlscf_.push_back(lnumber);

        getNumNonlinSolvIters(solver_memory_.get(), &lnumber);
        nnlsi_.push_back(lnumber);

        getNumJacEvals(solver_memory_.get(), &lnumber);
        nje_.push_back(lnumber);

        size_t number;
        getLastOrder(solver_memory_.get(), &number);
        order_.push_back(number);
    }

    void Solver::applyTolerances() const {
        if (!getInitDone())
            throw SunException("Solver instance was not yet set up, the "
                               "tolerances cannot be applied yet!");

        setSStolerances(rtol_, atol_);
    }

    void Solver::applyTolerancesFSA() const {
        if (!getInitDone())
            throw SunException("Solver instance was not yet set up, the "
                               "tolerances cannot be applied yet!");

        if (sensi_ < SensitivityOrder::first)
            return;

        if (np()) {
            std::vector<realtype> atols(np(), getAbsoluteToleranceFSA());
            setSensSStolerances(getRelativeToleranceFSA(), atols.data());
            setSensErrCon(true);
        }
    }

    double Solver::getRelativeToleranceFSA() const {
        return static_cast<double>(std::isnan(rtol_fsa_) ? rtol_ : rtol_fsa_);
    }

    void Solver::setRelativeToleranceFSA(const realtype rtol) {
        if (rtol < 0)
            throw SunException("rtol must be a non-negative number");

        rtol_fsa_ = rtol;

        if (getInitDone()) {
            applySensitivityTolerances();
        }
    }

    void Solver::initializeLinearSolver(const Model* model) const {

        (void) model;

        // Assert such that no default case is needed in switch statement.
        assert(linsol_ == LinearSolver::dense || linsol_ == LinearSolver::SuperLU);

        switch (linsol_) {

            case LinearSolver::dense:
                linear_solver_ = std::make_unique<SUNLinSolDense>(x_);
                setLinearSolver();
                setDenseJacFn();
                break;

            case LinearSolver::SuperLU:
                linear_solver_ = std::make_unique<SUNLinSolSuperLU>(
                        x_, model->nnz, CSC_MAT);
                setLinearSolver();
                setSparseJacFn();
                break;
        }
    }

    void Solver::initializeNonLinearSolver() const {

        assert(iter_== NonlinearSolverIteration::newton || iter_ == NonlinearSolverIteration::fixedpoint);

        switch (iter_) {
            case NonlinearSolverIteration::newton:
                non_linear_solver_ = std::make_unique<SUNNonLinSolNewton>(x_.getNVector());
                break;
            case NonlinearSolverIteration::fixedpoint:
                non_linear_solver_ =
                        std::make_unique<SUNNonLinSolFixedPoint>(x_.getNVector());
                break;
        }
        setNonLinearSolver();
    }

    bool Solver::getSensInitDone() const { return sens_initialized_; }

    void Solver::setSensInitDone() const { sens_initialized_ = true; }

    realtype Solver::getAbsoluteToleranceFSA() const {
        return std::isnan(atol_fsa_) ? atol_ : atol_fsa_;
    }

    void Solver::setAbsoluteToleranceFSA(const realtype atol) {
        if (atol < 0)
            throw SunException("atol must be a non-negative number");

        atol_fsa_ = atol;

        if (getInitDone()) {
            applySensitivityTolerances();
        }
    }

    void Solver::applySensitivityTolerances() const {
        if (sensi_ < SensitivityOrder::first)
            return;

        if (sensi_meth_ == SensitivityMethod::forward)
            applyTolerancesFSA();

        /*
        else if (sensi_meth_ == SensitivityMethod::adjoint && getAdjInitDone()) {
            for (int iMem = 0; iMem < (int)solver_memory_B_.size(); ++iMem)
                applyTolerancesASA(iMem);
        }
        */
    }

    SensitivityOrder Solver::getSensitivityOrder() const { return sensi_; }

    void Solver::setSensitivityOrder(const SensitivityOrder sensi) {
        if (sensi_ != sensi)
            (void) sensi_;
            // resetMutableMemory(nx(), nplist(), nquad());
        sensi_ = sensi;

        if (getInitDone())
            applySensitivityTolerances();
    }

    SensitivityMethod Solver::getSensitivityMethod() const { return sensi_meth_; }

    void Solver::setSensitivityMethod(const SensitivityMethod sensi_meth) {
        checkSensitivityMethod(sensi_meth);
        this->sensi_meth_ = sensi_meth;
    }

    void Solver::checkSensitivityMethod(const SensitivityMethod sensi_meth) const {
        if (sensi_meth != sensi_meth_)
            (void) sensi_meth;
            //resetMutableMemory(nx(), nplist(), nquad());
    }

    const Vector& Solver::getState(const realtype t) const {
        if (t == t_)
            return x_;

        if (solver_was_called_F_)
            getDky(t, 0);

        return dky_;
    }

    const Vector &Solver::getDerivativeState(const realtype t) const {
        if (t == t_)
            return dx_;

        if (solver_was_called_F_)
            getDky(t, 1);

        return dky_;
    }

    const VectorArray &Solver::getStateSensitivity(const realtype t) const {
        if (sens_initialized_ && solver_was_called_F_) {
            if (t == t_) {
                getSens();
            } else {
                getSensDky(t, 0);
            }
        }
        return sx_;
    }

    void wrapErrHandlerFn(int error_code, const char *module,
                          const char *function, char *msg, void* eh_data) {
        constexpr int BUF_SIZE = 250;
        char buffer[BUF_SIZE];
        char buffid[BUF_SIZE];
        snprintf(buffer, BUF_SIZE, "SUNEIGEN ERROR: in module %s in function %s : %s ", module,
                 function, msg);
        switch (error_code) {
            case 99:
                snprintf(buffid, BUF_SIZE, "SUNEIGEN:%s:%s:WARNING", module, function);
                break;

            case -1:
                snprintf(buffid, BUF_SIZE, "SUNEIGEN:%s:%s:TOO_MUCH_WORK", module, function);
                break;

            case -2:
                snprintf(buffid, BUF_SIZE, "SUNEIGEN:%s:%s:TOO_MUCH_ACC", module, function);
                break;

            case -3:
                snprintf(buffid, BUF_SIZE, "SUNEIGEN:%s:%s:ERR_FAILURE", module, function);
                break;

            case -4:
                snprintf(buffid, BUF_SIZE, "SUNEIGEN:%s:%s:CONV_FAILURE", module, function);
                break;

            default:
                snprintf(buffid, BUF_SIZE, "SUNEIGEN:%s:%s:OTHER", module, function);
                break;
        }


        if(!eh_data) {
            throw std::runtime_error("eh_data unset");
        }
        auto solver = static_cast<Solver const*>(eh_data);
        (void) solver;
        SunApplication app = SunApplication();
        app.warning(buffid, buffer);
    }

}  // namespace suneigen
