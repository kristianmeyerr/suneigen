#include "forward_problem.h"
#include "defines.h"
#include "solver.h"
#include "model.h"

namespace suneigen {

    ForwardProblem::ForwardProblem(Model* model_, Solver* solver_)
            : model(model_),
              solver(solver_),
              t_(model->t0()),
              x_(model->nx),
              x_old_(model->nx),
              dx_(model->nx),
              dx_old_(model->nx),
              xdot_(model->nx),
              xdot_old_(model->nx),
              sx_(model->nx,model->np()),
              sdx_(model->nx,model->np())
              {}

    void ForwardProblem::workForwardProblem() {
        FinalStateStorer fss(this);

        model->initialize(x_, dx_, sx_, sdx_,
                          solver->getSensitivityOrder() >=
                          SensitivityOrder::first);

        auto t0 = model->t0();
        solver->setup(t0, model, x_, dx_, sx_, sdx_);

        // update x0 after computing consistence IC/reinitialization
        x_ = solver->getState(model->t0());

        /* when computing forward sensitivities, we generally want to update sx
         after presimulation/preequilibration, and if we didn't do either this also
         wont harm.
        */
        if (solver->computingFSA() )
            sx_ = solver->getStateSensitivity(model->t0());

        /* store initial state and sensitivity*/
        initial_state_ = getSimulationState();

        /* loop over timepoints */
        for (it_ = 0; it_ < model->nt(); it_++) {
            auto nextTimepoint = model->getTimepoint(it_);

            if (std::isinf(nextTimepoint))
                break;

            if (nextTimepoint > model->t0()) {
                // Solve for nextTimepoint
                while (t_ < nextTimepoint) {
                    int status = solver->run(nextTimepoint);
                    solver->writeSolution(&t_, x_, dx_, sx_, dx_);
                    /* sx will be copied from solver on demand if sensitivities are computed */
                    (void)status;
                }
            }
            handleDataPoint(it_);
        }
    }

    void ForwardProblem::handleDataPoint(size_t /*it*/) {
        /* We only store the simulation state if it's not the initial state, as the
           initial state is stored anyway and we want to avoid storing it twice */
        if (t_ != model->t0() && timepoint_states_.count(t_) == 0)
            timepoint_states_[t_] = getSimulationState();
        /* store diagnosis information for debugging */
        solver->storeDiagnosis();
    }

    const SimulationState& ForwardProblem::getSimulationStateTimepoint(size_t it) const {
        if (model->getTimepoint(it) == initial_state_.t)
            return getInitialSimulationState();
        return timepoint_states_.find(model->getTimepoint(it))->second;
    }

    SimulationState ForwardProblem::getSimulationState() const {
        auto state = SimulationState();
        state.t = t_;
        state.x = x_;
        state.dx = dx_;
        if (solver->computingFSA() || t_ == model->t0())
            state.sx = sx_;

        return state;
    }

}
