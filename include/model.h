#ifndef SUNEIGEN_MODEL_H
#define SUNEIGEN_MODEL_H

#include "abstract_model.h"
#include "defines.h"
#include "sundials_matrix_wrapper.h"
#include "vector.h"
#include "model_dimensions.h"
#include "model_state.h"

namespace suneigen {

    class Model;
    class Solver;
    class SunApplication;

    class Model : public AbstractModel, public ModelDimensions {

    public:

        /** Default constructor */
        Model() = default;

        Model(ModelDimensions const& model_dimensions);

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
        using AbstractModel::fx0;

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
         * @brief Get total number of model parameters.
         * @return Length of parameter vector
         */
        [[nodiscard]] size_t np() const;

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
         * @brief Compute/get initial states.
         * @param x Output buffer.
         */
        void fx0(Vector &x);

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

        /**
         * Storage for model quantities beyond ModelState for the current timepoint
         */
        ModelStateDerived derived_state_;

        /** starting time */
        realtype tstart_ {0.0};

        SunApplication *app;

    protected:

        /** state initialization */
        std::vector<realtype> x0data_;

        /** sensitivity initialization (size nx_rdata x nplist, row-major) */
        std::vector<realtype> sx0data_;

        /** boolean indicating whether any entry in stateIsNonNegative is `true` */
        bool any_state_non_negative_ {false};

        Vector x_pos_tmp_ {0};

        /** vector of bools indicating whether state variables are to be assumed to
        * be positive */
        std::vector<bool> state_is_non_negative_;

        /**
         * @brief Timepoints for which model state/outputs/... are requested
         *
         * Vector of timepoints.
         */
        std::vector<realtype> ts_;

        /**
         * Indicates whether the result of every call to `Model::f*` should be
         * checked for finiteness
         */
        bool always_check_finite_ {false};

    };


}

#endif //SUNEIGEN_MODEL_H
