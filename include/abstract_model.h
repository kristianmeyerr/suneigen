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
         * @brief Root function
         * @param t time
         * @param x state
         * @param dx time derivative of state (DAE only)
         * @param root array to which values of the root function will be written
         */
        virtual void froot(realtype t, const Vector &x,
                           const Vector &dx, gsl::span<realtype> root) = 0;

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
         * @param p parameter vector
         * @param k constant vector
         */
        virtual void fx0(realtype *x0, realtype t, const realtype *p, const realtype *k);

        /**
         * @brief Model-specific implementation of fsx0
         * @param sx0 initial state sensitivities
         * @param t initial time
         * @param x0 initial state
         * @param ip sensitivity index
         */
        virtual void fsx0(realtype *sx0, realtype t, const realtype *x0,
                          const realtype *p, unsigned int ip);

        /**
         * @brief Model-specific implementation of fstau
         * @param stau total derivative of event timepoint
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param h Heaviside vector
         * @param sx current state sensitivity
         * @param ip sensitivity index
         * @param ie event index
         */
        virtual void fstau(realtype *stau, realtype t, const realtype *x,
                           const realtype *p, const realtype *h,
                           const realtype *sx, unsigned int ip, unsigned int ie);

        /**
         * @brief Model-specific implementation of fz
         * @param z value of event output
         * @param ie event index
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param k constant vector
         * @param h Heaviside vector
         */
        virtual void fz(realtype *z, unsigned int ie, realtype t, const realtype *x,
                        const realtype *p, const realtype *k, const realtype *h);

        /**
         * @brief Model-specific implementation of fdeltax
         * @param deltax state update
         * @param t current time
         * @param x current state
         * @param p parameter vector
         * @param h Heaviside vector
         * @param ie event index
         * @param xdot new model right hand side
         * @param xdot_old previous model right hand side
         */
        virtual void fdeltax(realtype *deltax, realtype t, const realtype *x,
                             const realtype *p, const realtype *h, unsigned int ie,
                             const realtype *xdot, const realtype *xdot_old);


    };

}

#endif //SUNEIGEN_ABSTRACT_MODEL_H
