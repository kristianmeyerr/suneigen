#ifndef SUNEIGEN_ABSTRACT_MODEL_H
#define SUNEIGEN_ABSTRACT_MODEL_H

#include "defines.h"
#include "vector.h"

#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include <memory>

namespace suneigen{

    class Solver;

    /**
     * @brief Abstract base class of suneigen::Model defining functions that need to
     * be implemented in a suneigen model.
     *
     * Some functions have empty default implementations or throw.
     * This class shall not have any data members.
     */
    class AbstractModel {
    public:
        virtual ~AbstractModel() = default;

        /**
         * @brief Retrieves the solver object
         * @return The Solver instance
         */
        virtual std::unique_ptr<Solver> getSolver() = 0;

        /**
         * @brief Residual function
         * @param t time
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @param xdot array to which values of the residual function will be
         * written
         */
        virtual void fxdot(realtype t, const Vector &x,
                           const Vector &dx, Vector &xdot) = 0;

        /**
         * @brief Dense Jacobian function
         * @param t time
         * @param cj scaling factor (inverse of timestep, DAE only)
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @param xdot values of residual function (unused)
         * @param J dense matrix to which values of the jacobian will be written
         */
        virtual void fJ(realtype t, realtype cj, const Vector &x,
                        const Vector &dx, const Vector &xdot,
                        SUNMatrix J) = 0;

        /**
         * @brief Sparse Jacobian function
         * @param t time
         * @param cj scaling factor (inverse of timestep, DAE only)
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @param xdot values of residual function (unused)
         * @param J sparse matrix to which values of the Jacobian will be written
         */
        virtual void fJSparse(realtype t, realtype cj, const Vector &x,
                              const Vector &dx, const Vector &xdot,
                              SUNMatrix J) = 0;

        /**
         * @brief Model-specific implementation of fx0
         * @param x0 initial state
         * @param t initial time
         */
        virtual void fx0(realtype *x0, realtype t) = 0;


    };

}

#endif //SUNEIGEN_ABSTRACT_MODEL_H
