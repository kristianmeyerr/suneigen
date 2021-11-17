#include "return_data.h"

#include "solver.h"

#include <cmath>

namespace suneigen {

    ReturnData::ReturnData(Solver const &solver, const Model &model)
            : ReturnData(model.getTimepoints(),
                         ModelDimensions(static_cast<ModelDimensions const&>(model)),
                         model.nt(), solver.getSensitivityOrder(),
                         solver.getSensitivityMethod())
                         {}

    ReturnData::ReturnData(std::vector<realtype> ts_,
                           ModelDimensions const& model_dimensions,
                           size_t nt_, SensitivityOrder sensi_,
                           SensitivityMethod sensi_meth_)
            : ModelDimensions(model_dimensions), ts(std::move(ts_)),
            nt(nt_), sensi(sensi_), sensi_meth(sensi_meth_),
            x_solver_(nx), sx_solver_(nx, np),
            x_rdata_(nx), sx_rdata_(nx, np){

        initiainitializeReporting();

    }

    void ReturnData::initiainitializeReporting(){

        xdot.resize(nx, std::nan(""));

        x.resize(nt * nx, 0.0);

        if (nt > 0) {
            numsteps.resize(nt, 0);
            numrhsevals.resize(nt, 0);
            numerrtestfails.resize(nt, 0);
            numnonlinsolvconvfails.resize(nt, 0);
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

        auto initialState = fwd.getInitialSimulationState();
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

        cpu_time = solver.getCpuTime();

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


}
