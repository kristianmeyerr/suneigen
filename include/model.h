#ifndef SUNEIGEN_MODEL_H
#define SUNEIGEN_MODEL_H

#include "abstract_model.h"
#include "defines.h"
#include "sundials_matrix_wrapper.h"
#include "simulation_parameters.h"
#include "vector.h"
#include "model_dimensions.h"
#include "model_state.h"

#include <map>
#include <memory>
#include <vector>

namespace suneigen {

    class Model;
    class Solver;

    /**
     * @brief The Model class represents a SunEigen ODE/DAE model.
     */
    class Model : public AbstractModel, public ModelDimensions {

    public:

        /** Default constructor */
        Model() = default;

        Model(ModelDimensions const& model_dimensions,
              SimulationParameters simulation_parameters,
              std::vector<int> z2event);

        /** Destructor. */
        ~Model() override = default;

        /**
         * @brief Copy assingment
         * @param other Model object.
         */
        Model(const Model& other) = delete;

        /**
         * @brief Copy assignment
         * @param other Object to copy from
         * @return
         */
        Model &operator=(Model const &other) = delete;

        // Overloaded base class methods
        using AbstractModel::fdeltax;
        using AbstractModel::fx0;
        using AbstractModel::fsx0;
        using AbstractModel::fstau;
        using AbstractModel::fz;

        /**
         * @brief Initialize model properties.
         * @param x Reference to state variables
         * @param dx Reference to time derivative of states (DAE only)
         * @param sx Reference to state variable sensitivities
         * @param sdx Reference to time derivative of state sensitivities (DAE only)
         * @param computeSensitivities Flag indicating whether sensitivities are to
         * be computed
         */
        void initialize(Vector &x, Vector &dx, VectorArray &sx,
                        VectorArray &sdx, bool computeSensitivities);

        /**
         * @brief Initialize initial states.
         * @param x State vector to be initialized
         */
        void initializeStates(Vector &x);

        /**
         * @brief Initialize initial state sensitivities.
         * @param sx Reference to state variable sensitivities
         * @param x Reference to state variables
         */
        void initializeStateSensitivities(VectorArray &sx, const Vector& x);

        /**
         * @brief Initialize the Heaviside variables `h` at the initial time `t0`.
         *
         * Heaviside variables activate/deactivate on event occurrences.
         *
         * @param x Reference to state variables
         * @param dx Reference to time derivative of states (DAE only)
         */
        void initHeaviside(const Vector &x, const Vector &dx);

        /**
         * @brief Get total number of model parameters.
         * @return Length of parameter vector
         */
        [[nodiscard]] size_t np() const;

        /**
         * @brief Get number of constants
         * @return Length of constant vector
         */
        [[nodiscard]] size_t nk() const;

        /**
        * @brief Get maximum number of events that may occur for each type.
        * @return Maximum number of events that may occur for each type
        */
        [[nodiscard]] unsigned int nMaxEvent() const;

        /**
         * @brief Get parameter vector.
         * @return The user-set parameters (see also `Model::getUnscaledParameters`)
         */
        [[nodiscard]] std::vector<realtype> const &getParameters() const;

        /**
         * @brief Set the parameter vector.
         * @param p Vector of parameters
         */
        void setParameters(std::vector<realtype> const &p);

        /**
         * @brief Set values for constants.
         * @param k Vector of fixed parameters
         */
        void setFixedParameters(std::vector<realtype> const &k);

        /**
         * @brief Set maximum number of events that may occur for each type.
         * @param nmaxevent Maximum number of events that may occur for each type
         */
        void setNMaxEvent(unsigned int nmaxevent);

        /**
         * @brief Get number of timepoints.
         * @return Number of timepoints
         */
        [[nodiscard]] size_t nt() const;

        /**
         * @brief Get the timepoint vector.
         * @return Timepoint vector
         */
        [[nodiscard]] std::vector<realtype> const& getTimepoints() const;

        /**
         * @brief Get simulation timepoint for time index `it`.
         * @param it Time index
         * @return Timepoint
         */
        [[nodiscard]] realtype getTimepoint(size_t it) const;

        /**
         * @brief Set the timepoint vector.
         * @param ts New timepoint vector
         */
        void setTimepoints(std::vector<realtype> const &ts);

        /**
         * @brief Get simulation start time.
         * @return Simulation start time
         */
        [[nodiscard]] double t0() const;

        /**
         * @brief Set simulation start time.
         * @param t0 Simulation start time
         */
        void setT0(double t0);

        /**
         * @brief Get flags indicating whether states should be treated as
         * non-negative.
         * @return Vector of flags
         */
        [[nodiscard]] std::vector<bool> const &getStateIsNonNegative() const;

        /**
         * @brief Set flags indicating whether states should be treated as
         * non-negative.
         * @param stateIsNonNegative Vector of flags
         */
        void setStateIsNonNegative(std::vector<bool> const &stateIsNonNegative);

        /**
         * @brief Set flags indicating that all states should be treated as
         * non-negative.
         */
        void setAllStatesNonNegative();

        /**
         * @brief Get the current model state.
         * @return Current model state
         */
        ModelState const &getModelState() const {
            return state_;
        }

        /**
         * @brief Set the current model state.
         * @param state Model state
         */
        void setModelState(ModelState const &state) {
            if (state.parameters.size() != np())
                throw SunException("Mismatch in parameter size");
            if (state.h.size() != ne)
                throw SunException("Mismatch in Heaviside size");
            if (state.x.size() != nx)
                throw SunException("Mismatch in conservation law size");
            state_ = state;
        }

        /**
         * @brief Get the initial states.
         * @return Initial state vector
         */
        std::vector<realtype> getInitialStates();

        /**
         * @brief Set the initial states.
         * @param x0 Initial state vector
         */
        void setInitialStates(std::vector<realtype> const &x0);

        /**
         * @brief Get the initial states sensitivities.
         * @return vector of initial state sensitivities
         */
        std::vector<realtype> getInitialStateSensitivities();

        /**
         * @brief Set the initial state sensitivities.
         * @param sx0 vector of initial state sensitivities with chainrule applied.
         * This could be a slice of ReturnData::sx or ReturnData::sx0
         */
        void setInitialStateSensitivities(std::vector<realtype> const &sx0);

        /**
         * @brief Get event-resolved observables.
         * @param z Output buffer (shape `nz`)
         * @param ie Event index
         * @param t Timepoint
         * @param x State variables
         */
        void getEvent(gsl::span<realtype> z, const unsigned int ie, const realtype t,
                      const Vector &x);

        /**
         * @brief Sensitivity of event timepoint, total derivative.
         *
         * Only forward sensitivities.
         *
         * @param stau Timepoint sensitivity (shape `nplist`)
         * @param t Timepoint
         * @param ie Event index
         * @param x State variables
         * @param sx State sensitivities
         */
        void getEventTimeSensitivity(std::vector<realtype> &stau, const realtype t,
                                     const unsigned int ie, const Vector &x,
                                     const VectorArray &sx);

        /**
         * @brief Update state variables after event.
         * @param x Current state (will be overwritten)
         * @param ie Event index
         * @param t Current timepoint
         * @param xdot Current residual function values
         * @param xdot_old Value of residual function before event
         */
        void addStateEventUpdate(Vector &x, const unsigned int ie, const realtype t,
                                 const Vector &xdot, const Vector &xdot_old);

        /**
         * @brief Update state sensitivity after event.
         * @param sx Current state sensitivity (will be overwritten)
         * @param ie Event index
         * @param t Current timepoint
         * @param x_old Current state
         * @param xdot Current residual function values
         * @param xdot_old Value of residual function before event
         * @param stau Timepoint sensitivity, to be computed with
         * `Model::getEventTimeSensitivity`
         */
        void addStateSensitivityEventUpdate(VectorArray &sx, const unsigned int ie,
                                            const realtype t,
                                            const Vector &x_old,
                                            const Vector &xdot,
                                            const Vector &xdot_old,
                                            const std::vector<realtype> &stau);

        /**
         * @brief Set the initial state sensitivities.
         * @param sx0 Vector of initial state sensitivities without chainrule
         * applied. This could be the readin from a `model.sx0data` saved to HDF5.
         */
        void setUnscaledInitialStateSensitivities(std::vector<realtype> const &sx0);

        /**
         * @brief Update the Heaviside variables `h` on event occurrences.
         *
         * @param rootsfound Provides the direction of the zero-crossing, so adding
         * it will give the right update to the Heaviside variables (zero if no root
         * was found)
         */
        void updateHeaviside(const std::vector<int> &rootsfound);

        /**
         * @brief Check if the given array has only finite elements.
         *
         * If not, try to give hints by which other fields this could be caused.
         *
         * @param array Array to check
         * @param fun Name of the function that generated the values (for more
         * informative messages).
         * @return `suneigen::SUNEIGEN_RECOVERABLE_ERROR` if a NaN/Inf value was found,
         * `amici::SUNEIGEN_SUCCESS` otherwise
         */
        int checkFinite(gsl::span<const realtype> array, const char* fun) const;

        /**
         * @brief Set whether the result of every call to `Model::f*` should be
         * checked for finiteness.
         * @param alwaysCheck bool
         */
        void setAlwaysCheckFinite(bool alwaysCheck);

        /**
         * @brief Get setting of whether the result of every call to `Model::f*`
         * should be checked for finiteness.
         * @return that
         */
        [[nodiscard]] bool getAlwaysCheckFinite() const;

        /**
         * @brief Compute/get initial states.
         * @param x Output buffer.
         */
        void fx0(Vector &x);

        /**
         * @brief Compute/get initial value for initial state sensitivities.
         * @param sx Output buffer for state sensitivities
         * @param x State variables
         */
        void fsx0(VectorArray &sx, const Vector &x);

        /**
         * @brief Copy states from solver to return_data
         * @param x_rdata Output buffer for state variables (stored in `suneigen::ReturnData`).
         * @param x_solver State variables (solver returns this)
         */
        void fx_rdata(Vector& x_rdata, const Vector& x_solver);

        /**
         * @brief Copy state sensitivites from solver to rdata
         * @param sx_rdata Output buffer for state variables sensitivities (stored in `suneigen::ReturnData`).
         * @param sx_solver State variables sensitivities (solver returns this)
         */
        void fsx_rdata(VectorArray& sx_rdata, const VectorArray& sx_solver);

        /**
         * @brief Compute fx_solver.
         *
         * @param x_solver State variables
         * @param x_rdata State variables
         */
        void fx_solver(realtype* x_solver, const realtype* x_rdata);

        /** Flag array for DAE equations */
        std::vector<realtype> idlist;

    protected:

        /**
         * @brief Write part of a slice to a buffer according to indices specified
                * in z2event.
        * @param slice Input data slice
        * @param buffer Output data slice
        * @param ie Event index
        */
        void writeSliceEvent(gsl::span<const realtype> slice,
                             gsl::span<realtype> buffer, int ie);

        /**
         * @brief Write part of a sensitivity slice to a buffer according to
         * indices specified in z2event.
         * @param slice source data slice
         * @param buffer output data slice
         * @param ie event index
         */
        void writeSensitivitySliceEvent(gsl::span<const realtype> slice,
                                        gsl::span<realtype> buffer, int ie);

        /**
         * @brief Compute event-resolved output.
         * @param ie Event index
         * @param t Current timepoint
         * @param x Current state
         */
        void fz(unsigned int ie, realtype t, const Vector &x);

        /**
         * @brief Compute fsx_solver.
         *
         * To be implemented by derived class if applicable.
         *
         * @param sx_rdata State sensitivity variables.
         * @param sx_solver State sensitivity variables.
         */
        virtual void fsx_solver(realtype *sx_solver, const realtype *sx_rdata);

        /**
         * @brief Compute non-negative state vector.
         *
         * Compute non-negative state vector according to stateIsNonNegative.
         * If anyStateNonNegative is set to `false`, i.e., all entries in
         * stateIsNonNegative are `false`, this function directly returns `x`,
         * otherwise all entries of x are copied in to `amici::Model::x_pos_tmp_`
         * and negative values are replaced by `0` where applicable.
         *
         * @param x State vector possibly containing negative values
         * @return State vector with negative values replaced by `0` according to
         * stateIsNonNegative
         */
        const_N_Vector computeX_pos(const_N_Vector x);

        /** All variables necessary for function evaluation */
        ModelState state_;

        /**
         * Storage for model quantities beyond ModelState for the current timepoint
         */
        ModelStateDerived derived_state_;

        /** index indicating to which event an event output belongs */
        std::vector<int> z2event_;

        /** state initialization */
        std::vector<realtype> x0data_;

        /** sensitivity initialization (size nx_rdata x nplist, row-major) */
        std::vector<realtype> sx0data_;

        /** vector of bools indicating whether state variables are to be assumed to
        * be positive */
        std::vector<bool> state_is_non_negative_;

        /** boolean indicating whether any entry in stateIsNonNegative is `true` */
        bool any_state_non_negative_ {false};

        /** maximal number of events to track */
        unsigned int nmaxevent_ {10};

        /**
         * Indicates whether the result of every call to `Model::f*` should be
         * checked for finiteness
         */
        bool always_check_finite_ {false};

    private:

        /** Simulation parameters, initial state, etc. */
        SimulationParameters simulation_parameters_;

    };


}

#endif //SUNEIGEN_MODEL_H
