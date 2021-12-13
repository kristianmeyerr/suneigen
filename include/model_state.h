#ifndef SUNEIGEN_MODEL_STATE_H
#define SUNEIGEN_MODEL_STATE_H

#include "defines.h"
#include "sundials_matrix_wrapper.h"
#include "model_dimensions.h"

#include <vector>

namespace suneigen {

    /**
     * @brief Exchange format to store and transfer the state of the
     * model at a specific timepoint.
     *
     * This is designed to only encompass the minimal
     * number of attributes that need to be transferred.
     */
    struct ModelState {
        /**
         * Flag indicating whether a certain Heaviside function should be active or
         * not (dimension: `ne`)
         */
        std::vector<realtype> h;

        /** States (dimension: `nx`) */
        std::vector<realtype> x;

        /** parameters (dimension: `np`) */
        std::vector<realtype> parameters;

        /** Constants (dimension: `nk`) */
        std::vector<realtype> fixedParameters;

    };

    /**
     * @brief Storage for `suneigen::Model` quantities computed for a specific timepoint.
     *
     * Serves as workspace for a model simulation to avoid repeated reallocation.
     */
    struct ModelStateDerived {
        ModelStateDerived() = default;

        /**
         * @brief Constructor from model dimensions.
         * @param dim Model dimensions
         */
        explicit ModelStateDerived(ModelDimensions const& dim);

        /** Sparse Jacobian (dimension: `amici::Model::nnz`) */
        SUNMatrixWrapper J_;

        /** temporary storage for flattened sx,
         * (dimension: `nx_solver` x `np`, row-major)
         */
        std::vector<realtype> sx_;

        /** temporary storage for event-resolved observable (dimension: nz) */
        std::vector<realtype> z_;

        /** temporary storage for `x_rdata` */
        std::vector<realtype> x_rdata_;

        /** temporary storage for `sx_rdata` slice */
        std::vector<realtype> sx_rdata_;

        /** temporary storage of positified state variables according to
         * stateIsNonNegative */
        Vector x_pos_tmp_ {0};

        /** temporary storage for change in x after event (dimension: `nx`) */
        std::vector<realtype> deltax_;
    };
}

#endif //SUNEIGEN_MODEL_STATE_H
