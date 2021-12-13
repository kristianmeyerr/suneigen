#ifndef SUNEIGEN_FORWARD_PROBLEM_H
#define SUNEIGEN_FORWARD_PROBLEM_H

#include "defines.h"
#include "vector.h"
#include "model.h"
#include "misc.h"

#include <map>

#include <iostream>

namespace suneigen {

    class Model;
    class Solver;
    class FinalStateStorer;

    /**
     * @brief implements an exchange format to store and transfer the state of a simulation at a
     * specific timepoint.
     */
    struct SimulationState{
        /** timepoint */
        realtype t;
        /** state variables */
        Vector x;
        /** state variables */
        Vector dx;
        /** state variable sensitivity */
        VectorArray sx;
        /** state of the model that was used for simulation */
        ModelState state;
    };

    /**
     * @brief The ForwardProblem class groups all functions for solving the
     * forward problem.
     */
    class ForwardProblem {
    public:
        /**
         * @brief Constructor
         * @param model pointer to Model instance
         * @param solver pointer to Solver instance
         * pass nullptr for no initialization
         */
        ForwardProblem(Model* model, Solver *solver);

        ~ForwardProblem() = default;

        /** allow FinalStateStorer to access private members and functions */
        friend ::suneigen::FinalStateStorer;

        /**
         * @brief Solve the forward problem.
         *
         * If forward sensitivities are enabled this will also compute sensitivities.
         */
        void workForwardProblem();

        /**
         * @brief Accessor for t
         * @return t
         */
        [[nodiscard]] realtype getTime() const {
            return t_;
        }

        /**
         * @brief Accessor for x
         * @return x
         */
        [[nodiscard]] Vector const& getState() const {
            return x_;
        }

        /**
         * @brief Accessor for dx
         * @return dx
         */
        [[nodiscard]] Vector const& getStateDerivative() const {
            return dx_;
        }

        /**
         * @brief Accessor for sx
         * @return sx
         */
        [[nodiscard]] VectorArray const& getStateSensitivity() const {
            return sx_;
        }

        /**
         * @brief Accessor for x_disc
         * @return x_disc
         */
        [[nodiscard]] std::vector<Vector> const& getStatesAtDiscontinuities() const {
            return x_disc_;
        }

        /**
         * @brief Accessor for xdot_disc
         * @return xdot_disc
         */
        [[nodiscard]] std::vector<Vector> const& getRHSAtDiscontinuities() const {
            return xdot_disc_;
        }

        /**
         * @brief Accessor for xdot_old_disc
         * @return xdot_old_disc
         */
        [[nodiscard]] std::vector<Vector> const& getRHSBeforeDiscontinuities() const {
            return xdot_old_disc_;
        }

        /**
         * @brief Accessor for nroots
         * @return nroots
         */
        [[nodiscard]] std::vector<unsigned int> const& getNumberOfRoots() const {
            return nroots_;
        }

        /**
         * @brief Accessor for discs
         * @return discs
         */
        [[nodiscard]] std::vector<realtype> const& getDiscontinuities() const {
            return discs_;
        }

        /**
         * @brief Accessor for rootidx
         * @return rootidx
         */
        [[nodiscard]] std::vector<std::vector<int>> const& getRootIndexes() const {
            return root_idx_;
        }

        /**
         * @brief Accessor for pointer to x
         * @return &x
         */
        Vector* getStatePointer() {
            return &x_;
        }

        /**
         * @brief Accessor for pointer to dx
         * @return &dx
         */
        Vector *getStateDerivativePointer() {
            return &dx_;
        }

        /**
         * @brief accessor for pointer to sx
         * @return &sx
         */
        VectorArray *getStateSensitivityPointer() {
            return &sx_;
        }

        /**
         * @brief Accessor for pointer to sdx
         * @return &sdx
         */
        VectorArray *getStateDerivativeSensitivityPointer() {
            return &sdx_;
        }

        /**
         * @brief Accessor for it
         * @return it
         */
        [[nodiscard]] size_t getCurrentTimeIteration() const {
            return it_;
        }

        /**
         * @brief Returns final time point for which simulations are available
         * @return time point
         */
        [[nodiscard]] realtype getFinalTime() const {
            return final_state_.t;
        }

        /**
         * @brief Returns maximal event index for which simulations are available
         * @return index
         */
        [[nodiscard]] int getEventCounter() const {
            return static_cast<int>(event_states_.size() - 1);
        }

        /**
         * @brief Returns maximal event index for which the timepoint is available
         * @return index
         */
        [[nodiscard]] unsigned int getRootCounter() const {
            return static_cast<unsigned int>(discs_.size() - 1);
        }

        /**
         * @brief Retrieves the carbon copy of the simulation state variables at
         * the specified timepoint index
         * @param it timepoint index
         * @return state
         */
        [[nodiscard]] const SimulationState& getSimulationStateTimepoint(size_t it) const;

        /**
         * @brief Retrieves the carbon copy of the simulation state variables at
         * the specified event index
         * @param iroot event index
         * @return SimulationState
         */
        [[nodiscard]] const SimulationState &getSimulationStateEvent(unsigned int iroot) const {
            return event_states_.at(iroot);
        }

        /**
         * @brief Retrieves the carbon copy of the simulation state variables at the
         * initial timepoint
         * @return SimulationState
         */
        [[nodiscard]] const SimulationState& getInitialSimulationState() const {
            return initial_state_;
        }

        /**
         * @brief Retrieves the carbon copy of the simulation state variables at the
         * final timepoint (or when simulation failed)
         * @return SimulationState
         */
        [[nodiscard]] const SimulationState& getFinalSimulationState() const {
            return final_state_;
        }

        /** pointer to model instance */
        Model* model;

        /** pointer to solver instance */
        Solver* solver;

    private:

        /**
         * @brief Execute everything necessary for the handling of events
         *
         * @param tlastroot pointer to the timepoint of the last event
         * @param seflag Secondary event flag
         */
        void handleEvent(realtype *tlastroot, bool seflag);

        /**
         * @brief Extract output information for events
         */
        void storeEvent();

        /**
         * @brief Execute everything necessary for the handling of data points
         *
         * @param it index of data point
         */
        void handleDataPoint(size_t it);

        /**
         * @brief Applies the event bolus to the current state
         */
        void applyEventBolus();

        /**
         * @brief Applies the event bolus to the current sensitivities
         */
        void applyEventSensiBolusFSA();

        /**
         * @brief checks whether there are any events to fill
         *
         * @param nmaxevent maximal number of events
         */
        [[nodiscard]] bool checkEventsToFill(int nmaxevent) const {
            return std::any_of(nroots_.cbegin(), nroots_.cend(),
                               [nmaxevent](int curNRoots) {
                                   return curNRoots < nmaxevent;
                               });
        }

        /**
         * @brief fills events at final timepoint if necessary
         *
         * @param nmaxevent maximal number of events
         */
        void fillEvents(int nmaxevent) {
            if (checkEventsToFill(nmaxevent)) {
                discs_.push_back(t_);
                storeEvent();
            }
        }

        /**
         * @brief Creates a carbon copy of the current simulation state variables
         * @return state
         */
        [[nodiscard]] SimulationState getSimulationState() const;

        /** array of index vectors (dimension ne) indicating whether the respective
         * root has been detected for all so far encountered discontinuities,
         * extended as needed (dimension: dynamic) */
        std::vector<std::vector<int>> root_idx_;

        /** array of number of found roots for a certain event type
         * (dimension: ne) */
        std::vector<unsigned int> nroots_;

        /** array of values of the root function (dimension: ne) */
        std::vector<realtype> rootvals_;

        /** temporary rootval storage to check crossing in secondary event
        * (dimension: ne) */
        std::vector<realtype> rval_tmp_;

        /** array containing the time-points of discontinuities
        * (dimension: nmaxevent x ne, ordering = ?) */
        std::vector<realtype> discs_;

        /** array containing the index of discontinuities
        * (dimension: nmaxevent x ne, ordering = ?) */
        std::vector<realtype> irdiscs_;

        /** array of state vectors (dimension nx) for all so far encountered
         * discontinuities, extended as needed (dimension dynamic) */
        std::vector<Vector> x_disc_;

        /** array of differential state vectors (dimension nx) for all so far
        *   encountered discontinuities, extended as needed (dimension dynamic) */
        std::vector<Vector> xdot_disc_;

        /** array of old differential state vectors (dimension nx) for all so far
         * encountered discontinuities, extended as needed (dimension dynamic) */
        std::vector<Vector> xdot_old_disc_;

        /** current time */
        realtype t_;

        /**
         * @brief Array of flags indicating which root has been found.
         *
         * Array of length nr (ne) with the indices of the user functions gi found
         * to have a root. For i = 0, . . . ,nr 1 if gi has a root, and = 0 if not.
         */
        std::vector<int> roots_found_;

        /** simulation states history at timepoints  */
        std::map<realtype, SimulationState> timepoint_states_;

        /** simulation state history at events*/
        std::vector<SimulationState> event_states_;

        /** simulation state after initialization*/
        SimulationState initial_state_;

        /** simulation state after simulation*/
        SimulationState final_state_;

        /** state vector */
        Vector x_;

        /** old state vector */
        Vector x_old_;

        /** differential state vector */
        Vector dx_;

        /** old differential state vector */
        Vector dx_old_;

        /** time derivative state vector */
        Vector xdot_;

        /** old time derivative state vector */
        Vector xdot_old_;

        /** sensitivity state vector array */
        VectorArray sx_;

        /** differential sensitivity state vector array */
        VectorArray sdx_;

        /** sensitivity of the event timepoint (dimension: np) */
        std::vector<realtype> stau_;

        /** storage for last found root */
        realtype tlastroot_ {0.0};

        /** current iteration number for time index */
        size_t it_{0};

    };

    /**
     * @brief stores the stimulation state when it goes out of scope
     */
    class FinalStateStorer : public ContextManager {
    public:
        /**
         * @brief constructor, attaches problem pointer
         * @param fwd problem from which the simulation state is to be stored
         */
        explicit FinalStateStorer(ForwardProblem* fwd) : fwd_(fwd) {}

        FinalStateStorer &operator=(const FinalStateStorer &other) = delete;

        /**
         * @brief destructor, stores simulation state
         */
        ~FinalStateStorer() {
            if(fwd_)
                fwd_->final_state_ = fwd_->getSimulationState();
        }
    private:
        ForwardProblem* fwd_;
    };

}

#endif //SUNEIGEN_FORWARD_PROBLEM_H
