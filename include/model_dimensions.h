#ifndef SUNEIGEN_MODEL_DIMENSIONS_H
#define SUNEIGEN_MODEL_DIMENSIONS_H

#include <gsl/gsl-lite.hpp>
#include <vector>

namespace suneigen {

    /**
     * @brief Container for model dimensions.
     *
     * Holds number of states, observables, etc.
     */
    struct ModelDimensions {
        /** Default ctor */
        ModelDimensions() = default;

        ModelDimensions(
                const size_t nx_,
                const size_t np_,
                const size_t nnz_)
                : nx(nx_), np(np_), nnz(nnz_) {}

        /** Number of states */
        size_t nx{0};

        /** Number of parameters */
        size_t np{0};

        /** Number of nonzero entries in Jacobian */
        size_t nnz{0};
    };
}

#endif //SUNEIGEN_MODEL_DIMENSIONS_H
