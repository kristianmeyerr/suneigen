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

        ModelDimensions() = default;

        /**
         * @brief Constructor with model dimensions
         * @param nx_ Number of state variables
         * @param np_ Number of parameters
         * @param nk_ Number of constants
         * @param nz_ Number of event observables
         * @param ny_ Number of observables
         * @param ne_ Number of events
         * @param nnz_ Number of nonzero elements in Jacobian
         * @param ubw_ Upper matrix bandwidth in the Jacobian
         * @param lbw_ Lower matrix bandwidth in the Jacobian
         */
        ModelDimensions(
                const size_t nx_,
                const size_t np_,
                const size_t nk_,
                const size_t ny_,
                const size_t nz_,
                const unsigned int ne_,
                const size_t nnz_,
                const size_t ubw_,
                const size_t lbw_)
                : nx(nx_), np(np_), nk(nk_), ny(ny_), nz(nz_), ne(ne_), nnz(nnz_), ubw(ubw_), lbw(lbw_)
                {
                }

        /** Number of states */
        size_t nx{0};

        /** Number of parameters */
        size_t np{0};

        /** Number of constants */
        size_t nk{0};

        /** Number of observables */
        size_t ny{0};

        /** Number of event outputs */
        size_t nz{0};

        /** Number of events */
        unsigned int ne{0};

        /** Number of nonzero entries in Jacobian */
        size_t nnz{0};

        /** Upper bandwidth of the Jacobian */
        size_t ubw{0};

        /** Lower bandwidth of the Jacobian */
        size_t lbw{0};

    };
}

#endif //SUNEIGEN_MODEL_DIMENSIONS_H
