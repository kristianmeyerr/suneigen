#include "return_data.h"
#include "solver.h"
#include "defines.h"

#include <string>
#include <iostream>
#include <cmath>

namespace suneigen {

    ReturnData::ReturnData(Solver const &solver, const Model &model)
            : ReturnData(model.getTimepoints(),
                         ModelDimensions(static_cast<ModelDimensions const&>(model)),
                         model.nMaxEvent(), model.nt(), solver.getSensitivityOrder(),
                         solver.getSensitivityMethod())
                         {}

    ReturnData::ReturnData(std::vector<realtype> ts_,
                           ModelDimensions const& model_dimensions,
                           unsigned int nmaxevent_, size_t nt_, SensitivityOrder sensi_,
                           SensitivityMethod sensi_meth_)
            : ModelDimensions(model_dimensions), ts(std::move(ts_)),
              nmaxevent(nmaxevent_), nt(nt_), sensi(sensi_), sensi_meth(sensi_meth_),
            x_solver_(nx), sx_solver_(nx, np),
            x_rdata_(nx), sx_rdata_(nx, np), nroots_(ne){

        initiainitializeReporting();

    }

    void ReturnData::initiainitializeReporting(){

        xdot.resize(nx, std::nan(""));

        // initialize with 0.0, so we only need to write non-zero values
        z.resize(nmaxevent * nz, 0.0);

        x.resize(nt * nx, 0.0);

        if (nt > 0) {
            numsteps.resize(nt, 0);
            numrhsevals.resize(nt, 0);
            numerrtestfails.resize(nt, 0);
            numnonlinsolvconvfails.resize(nt, 0);
            numnonlinsolviter.resize(nt, 0);
            numjacevals.resize(nt, 0);
            order.resize(nt, 0);
        }

        x0.resize(nx, std::nan(""));
        if (sensi >= SensitivityOrder::first) {
            sx0.resize(nx * np, std::nan(""));

            if (sensi_meth == SensitivityMethod::forward ||
                sensi >= SensitivityOrder::second) {
                // for second order we can fill in from the augmented states
                sx.resize(nt * nx * np, 0.0);
            }
        }

    }

    void ReturnData::processSimulationObjects(ForwardProblem const *fwd,
                                              Model &model, Solver const &solver){
        processSolver(solver);

        if (fwd)
            processForwardProblem(*fwd, model);
        else
            invalidate(0);
    }

    void ReturnData::processForwardProblem(ForwardProblem const &fwd, Model &model) {

        t0 = model.t0();
        SimulationState initialState = fwd.getInitialSimulationState();
        if (initialState.x.getLength() == 0 && model.nx > 0)
            return; // if x wasn't set forward problem failed during initialization

        readSimulationState(initialState, model);

        if (!x0.empty()) {
            model.fx_rdata(x_rdata_, x_solver_);
            writeSlice(x_rdata_, x0);
        }

        // process time point data
        realtype tf = fwd.getFinalTime();
        for (unsigned int it = 0; it < model.nt(); it++) {
            if (model.getTimepoint(it) <= tf) {
                readSimulationState(fwd.getSimulationStateTimepoint(it), model);
                getDataOutput(it, model);
            } else {
                // check for integration failure but consider postequilibration
                if (!std::isinf(model.getTimepoint(it)))
                    invalidate(it);
            }
        }

        // process event data
        if (nz > 0) {
            auto rootidx = fwd.getRootIndexes();
            for (int iroot = 0; iroot <= fwd.getEventCounter(); iroot++) {
                readSimulationState(fwd.getSimulationStateEvent(static_cast<unsigned int>(iroot)), model);
                getEventOutput(t_, rootidx.at(static_cast<unsigned int>(iroot)), model);
            }
        }

    }

    void ReturnData::getEventOutput(realtype t, std::vector<int> rootidx,
                                    Model &model) {

        for (unsigned int ie = 0; ie < ne; ie++) {
            if (rootidx.at(ie) != 1 || nroots_.at(ie) >= nmaxevent)
                continue;

            /* get event output */
            if (!z.empty())
                model.getEvent(slice(z, nroots_.at(ie), nz), ie, t, x_solver_);

            if (sensi >= SensitivityOrder::first) {
                if (sensi_meth == SensitivityMethod::forward) {
                    // getEventSensisFSA(ie, t, model);
                }
            }
            nroots_.at(ie)++;
        }
    }

    void ReturnData::invalidate(const size_t it_start) {
        if (it_start >= nt)
            return;

        if (!x.empty())
            std::fill(x.begin() + static_cast<long>(nx * it_start), x.end(), std::nan(""));

        if (!sx.empty())
            std::fill(sx.begin() + static_cast<long>(nx * np * it_start), sx.end(), std::nan(""));
    }

    void ReturnData::readSimulationState(SimulationState const& state,
                                         Model& model) {
        x_solver_ = state.x;
        dx_solver_ = state.dx;
        if (computingFSA() || state.t == model.t0())
            sx_solver_ = state.sx;
        t_ = state.t;
    }

    void ReturnData::getDataOutput(size_t it, Model &model) {

        if (!x.empty()) {
            model.fx_rdata(x_rdata_, x_solver_);
            writeSlice(x_rdata_, slice(x, it, nx));
        }

        if (sensi >= SensitivityOrder::first && np > 0) {

            if (sensi_meth == SensitivityMethod::forward)
                getDataSensisFSA(it, model);

        }
    }

    void ReturnData::processSolver(Solver const &solver) {

        // Solver statistics
        cpu_time = solver.getCpuTime();
        solver_type = solver.getSolverType();
        abstol = solver.getAbsoluteTolerance();
        reltol = solver.getRelativeTolerance();
        lmm_solver = solver.getLinearMultistepMethod();
        ls_solver = solver.getLinearSolver();

        const std::vector<size_t> *tmp;

        if (!numsteps.empty()) {
            tmp = &solver.getNumSteps();
            // copy_n instead of assignment to ensure length `nt`
            // (vector from solver may be shorter in case of integration errors)
            std::copy_n(tmp->cbegin(), tmp->size(), numsteps.begin());
        }

        if (!numsteps.empty()) {
            tmp = &solver.getNumRhsEvals();
            std::copy_n(tmp->cbegin(), tmp->size(), numrhsevals.begin());
        }

        if (!numerrtestfails.empty()) {
            tmp = &solver.getNumErrTestFails();
            std::copy_n(tmp->cbegin(), tmp->size(), numerrtestfails.begin());
        }

        if (!numnonlinsolvconvfails.empty()) {
            tmp = &solver.getNumNonlinSolvConvFails();
            std::copy_n(tmp->cbegin(), tmp->size(), numnonlinsolvconvfails.begin());
        }

        if (!numnonlinsolviter.empty()) {
            tmp = &solver.getNumNonlinSolvIters();
            std::copy_n(tmp->cbegin(), tmp->size(), numnonlinsolviter.begin());
        }

        if (!numjacevals.empty()) {
            tmp = &solver.getNumJacEvals();
            std::copy_n(tmp->cbegin(), tmp->size(), numjacevals.begin());
        }

        if (!order.empty()) {
            tmp = &solver.getLastOrder();
            std::copy_n(tmp->cbegin(), tmp->size(), order.begin());
        }
    }

    void ReturnData::getDataSensisFSA(size_t it, Model &model) {
        if (!sx.empty()) {
            model.fsx_rdata(sx_rdata_, sx_solver_);
            for (unsigned int ip = 0; ip < np; ip++) {
                writeSlice(sx_rdata_[ip],
                           slice(sx, it * np + ip, nx));
            }
        }
    }

    void ReturnData::printStatistics(){

        // Statistics

        int numsteps_total = 0;
        for (auto& n : numsteps)
            numsteps_total += n;

        int numrhsevals_total = 0;
        for (auto& n : numrhsevals)
            numrhsevals_total += n;

        int numerrtestfails_total = 0;
        for (auto& n : numerrtestfails)
            numerrtestfails_total += n;

        int numnonlinsolvconvfails_total = 0;
        for (auto& n : numnonlinsolvconvfails)
            numnonlinsolvconvfails_total += n;

        int numnonlinsolviter_total = 0;
        for (auto& n : numnonlinsolviter)
            numnonlinsolviter_total += n;

        int numjacevals_total = 0;
        for (auto& n : numjacevals)
            numjacevals_total += n;

        std::string lmm;
        switch(lmm_solver){
            case LinearMultistepMethod::BDF:
                lmm = "BDF";
                break;
            case LinearMultistepMethod::adams:
                lmm = "adams";
                break;
        }

        std::string ls;
        switch(ls_solver){
            case LinearSolver::dense:
                ls = "Dense";
                break;
            case LinearSolver::SuperLU:
                ls = "SuperLU";
                break;
        }

        int max_order = *max_element(std::begin(order), std::end(order));

        std::cout << "[[Simulation statistics]]" << std::endl;
        std::cout << "    # Number of steps: " << numsteps_total << std::endl;
        std::cout << "    # Number of function evaluations: " << numrhsevals_total << std::endl;
        std::cout << "    # Number of Jacobian evaluations: " << numjacevals_total << std::endl;
        std::cout << "    # Number of error test failures: " << numerrtestfails_total << std::endl;
        std::cout << "    # Number of nonlinear iterations: " << numnonlinsolviter_total << std::endl;
        std::cout << "    # Number of nonlinear convergence failures: " << numnonlinsolvconvfails_total << std::endl;
        std::cout << "    # Simulation interval: " << ts.back() - t0 << std::endl;
        std::cout << "    # Wall clock time (ms): " << cpu_time << std::endl;
        std::cout << std::endl;
        std::cout << "[[Solver options]]" << std::endl;
        std::cout << "    # Solver: " << solver_type << std::endl;
        std::cout << "    # Linear multistep method: " << lmm << std::endl;
        std::cout << "    # Linear solver: " << ls << std::endl;
        std::cout << "    # Maximal order: " << max_order << std::endl;
        std::cout << "    # Absolute tolerance: " << abstol << std::endl;
        std::cout << "    # Relative tolerance: " << reltol << std::endl;
        std::cout << std::endl;




    }


}
