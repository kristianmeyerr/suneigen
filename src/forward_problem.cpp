#include "forward_problem.h"
#include "defines.h"
#include "solver.h"
#include "model.h"
#include "suneigen.h"

#include <iostream>

namespace suneigen {

    ForwardProblem::ForwardProblem(Model* model_, Solver* solver_)
            : model(model_),
              solver(solver_),
              nroots_(model->ne, 0),
              rootvals_(model->ne, 0.0),
              rval_tmp_(model->ne, 0.0),
              t_(model->t0()),
              roots_found_(model->ne, 0),
              x_(model->nx),
              x_old_(model->nx),
              dx_(model->nx),
              dx_old_(model->nx),
              xdot_(model->nx),
              xdot_old_(model->nx),
              sx_(model->nx,model->np()),
              sdx_(model->nx,model->np()),
              stau_(model->np())
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
                    /* sx will be copied from solver on demand if sensitivities
                    are computed */
                    if (status == SUNEIGEN_ILL_INPUT) {
                        /* clustering of roots => turn off rootfinding */
                        solver->turnOffRootFinding();
                    } else if (status == SUNEIGEN_ROOT_RETURN) {
                        handleEvent(&tlastroot_, false);
                    }
                }
            }
            handleDataPoint(it_);
        }
    }

    void ForwardProblem::storeEvent() {

        if (t_ == model->getTimepoint(model->nt() - 1)) {
            // call from fillEvent at last timepoint
            model->froot(t_, x_, dx_, rootvals_);
            for (unsigned int ie = 0; ie < model->ne; ie++) {
                roots_found_.at(ie) = (nroots_.at(ie) < model->nMaxEvent()) ? 1 : 0;
            }
            root_idx_.push_back(roots_found_);
        }

        if (static_cast<int>(getRootCounter()) < getEventCounter()) {
            /* update stored state (sensi) */
            std::cout << getRootCounter() << std::endl;
            std::cout << getEventCounter() << std::endl;
            event_states_.at(getRootCounter()) = getSimulationState();
        } else {
            /* add stored state (sensi) */
            event_states_.push_back(getSimulationState());
        }

        /* EVENT OUTPUT */
        for (unsigned int ie = 0; ie < model->ne; ie++) {
            /* only look for roots of the rootfunction not discontinuities */
            if (nroots_.at(ie) >= model->nMaxEvent())
                continue;

            /* only consider transitions false -> true or event filling */
            if (roots_found_.at(ie) != 1 &&
                t_ != model->getTimepoint(model->nt() - 1)) {
                continue;
            }

            nroots_.at(ie)++;
        }

        if (t_ == model->getTimepoint(model->nt() - 1)) {
            // call from fillEvent at last timepoint
            // loop until all events are filled
            fillEvents(static_cast<int>(model->nMaxEvent()));
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

    void ForwardProblem::applyEventBolus() {
        for (unsigned int ie = 0; ie < model->ne; ie++)
            if (roots_found_.at(ie) == 1) // only consider transitions false -> true
                model->addStateEventUpdate(x_, ie, t_, xdot_, xdot_old_);
    }

    void ForwardProblem::applyEventSensiBolusFSA() {
        for (unsigned int ie = 0; ie < model->ne; ie++)
            if (roots_found_.at(ie) == 1) // only consider transitions false -> true
                /*  */
                model->addStateSensitivityEventUpdate(sx_, ie, t_, x_old_, xdot_,
                                                      xdot_old_, stau_);
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

    void ForwardProblem::handleEvent(realtype *tlastroot, const bool seflag) {

        /* store Heaviside information at event occurrence */
        model->froot(t_, x_, dx_, rootvals_);

        /* store timepoint at which the event occurred*/
        discs_.push_back(t_);

        /* extract and store which events occurred */
        if (!seflag) {
            solver->getRootInfo(roots_found_.data());
        }
        root_idx_.push_back(roots_found_);

        rval_tmp_ = rootvals_;

        if (!seflag) {
            /* only check this in the first event fired, otherwise this will always
             * be true */
            if (t_ == *tlastroot) {
                throw SunException("SunEigen is stuck in an event, as the initial "
                                   "step-size after the event is too small. "
                                   "To fix this, increase absolute and relative "
                                   "tolerances!");
            }
            *tlastroot = t_;
        }

        if(model->nz>0)
            storeEvent();

        /* if we need to do forward sensitivities later on we need to store the old
         * x and the old xdot */
        if (solver->getSensitivityOrder() >= SensitivityOrder::first) {
            /* store x and xdot to compute jump in sensitivities */
            x_old_ = x_;
        }
        if (solver->computingFSA()) {
            model->fxdot(t_, x_, dx_, xdot_);
            xdot_old_ = xdot_;
            dx_old_ = dx_;
            /* compute event-time derivative only for primary events, we get
             * into trouble with multiple simultaneously firing events here (but
             * is this really well defined then?), in that case just use the
             * last ie and hope for the best. */
            if (!seflag) {
                for (unsigned int ie = 0; ie < model->ne; ie++) {
                    if (roots_found_.at(ie) == 1) {
                        /* only consider transitions false -> true */
                        model->getEventTimeSensitivity(stau_, t_, ie, x_, sx_);
                    }
                }
            }
        } else if (solver->computingASA()) {
            /* store x to compute jump in discontinuity */
            x_disc_.push_back(x_);
            xdot_disc_.push_back(xdot_);
            xdot_old_disc_.push_back(xdot_old_);
        }

        model->updateHeaviside(roots_found_);

        applyEventBolus();

        if (solver->computingFSA()) {
            /* compute the new xdot  */
            model->fxdot(t_, x_, dx_, xdot_);
            applyEventSensiBolusFSA();
        }

        int secondevent = 0;

        /* check whether we need to fire a secondary event */
        model->froot(t_, x_, dx_, rootvals_);
        for (unsigned int ie = 0; ie < model->ne; ie++) {
            /* the same event should not trigger itself */
            if (roots_found_.at(ie) == 0) {
                /* check whether there was a zero-crossing */
                if (0 > rval_tmp_.at(ie) * rootvals_.at(ie)) {
                    if (rval_tmp_.at(ie) < rootvals_.at(ie)) {
                        roots_found_.at(ie) = 1;
                    } else {
                        roots_found_.at(ie) = -1;
                    }
                    secondevent++;
                } else {
                    roots_found_.at(ie) = 0;
                }
            } else {
                /* don't fire the same event again */
                roots_found_.at(ie) = 0;
            }
        }
        /* fire the secondary event */
        if (secondevent > 0) {
            /* Secondary events may result in wrong forward sensitivities,
             * if the secondary event has a bolus... */
            if (solver->computingFSA()){
                SunApplication app;
                app.warning("SUNEIGEN:simulation",
                            "Secondary event was triggered. Depending on "
                            "the bolus of the secondary event, forward "
                            "sensitivities can be incorrect.");
            }
            handleEvent(tlastroot, true);
        }

        /* only reinitialise in the first event fired */
        if (!seflag) {
            solver->reInit(t_, x_, dx_);
            if (solver->computingFSA()) {
                solver->sensReInit(sx_, sdx_);
            }
        }
    }

}
