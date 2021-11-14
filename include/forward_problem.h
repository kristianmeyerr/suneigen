#ifndef SUNEIGEN_FORWARD_PROBLEM_H
#define SUNEIGEN_FORWARD_PROBLEM_H

#include "defines.h"
#include "vector.h"
#include "misc.h"

#include <map>

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
         * @brief Retrieves the carbon copy of the simulation state variables at
         * the specified timepoint index
         * @param it timepoint index
         * @return state
         */
        const SimulationState& getSimulationStateTimepoint(size_t it) const;

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
         * @brief Creates a carbon copy of the current simulation state variables
         * @return state
         */
        [[nodiscard]] SimulationState getSimulationState() const;

        /**
         * @brief Execute everything necessary for the handling of data points
         *
         * @param it index of data point
         */
        void handleDataPoint(size_t it);

        /** current time */
        realtype t_;

        /** simulation states history at timepoints  */
        std::map<realtype, SimulationState> timepoint_states_;

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

        /** simulation state after initialization*/
        SimulationState initial_state_;

        /** simulation state after simulation*/
        SimulationState final_state_;

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
        ForwardProblem * fwd_;
    };

}

#endif //SUNEIGEN_FORWARD_PROBLEM_H
