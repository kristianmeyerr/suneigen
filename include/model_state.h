#ifndef SUNEIGEN_MODEL_STATE_H
#define SUNEIGEN_MODEL_STATE_H

#include "defines.h"
#include "sundials_matrix_wrapper.h"
#include "model_dimensions.h"

#include <vector>

namespace suneigen {

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

        /** temporary storage for `x_rdata` */
        std::vector<realtype> x_rdata_;

        /** temporary storage for `sx_rdata` slice */
        std::vector<realtype> sx_rdata_;

        /** temporary storage of positified state variables according to
         * stateIsNonNegative */
        Vector x_pos_tmp_ {0};
    };
}

#endif //SUNEIGEN_MODEL_STATE_H
