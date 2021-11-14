#ifndef SUNEIGEN_MODEL_ODE_H
#define SUNEIGEN_MODEL_ODE_H

#include "model.h"

namespace suneigen {

    class CVodeSolver;

    /**
     * @brief The Model class represents an SUNEIGEN ODE model.
     *
     * The model does not contain any data, but represents the state
     * of the model at a specific time t. The states must not always be
     * in sync, but may be updated asynchronously.
     */
    class Model_ODE : public Model {
    public:
        /** default constructor */
        Model_ODE() = default;

        /**
         * @brief Constructor with model dimensions
         * @param model_dimensions Model dimensions
         */
        explicit Model_ODE(ModelDimensions const& model_dimensions)
                : Model(model_dimensions) {}

        /**
         * @brief Implementation of fxdot using Vectors
         * @param t timepoint
         * @param x Vector with the states
         * @param xdot Vector with the right hand side
         */
        void fxdot(realtype t, const Vector &x, const Vector &dx,
                   Vector &xdot) override;

        /**
         * @brief Implementation of fxdot at the N_Vector level, this function
         * provides an interface to the model specific routines for the solver
         * implementation as well as the Vector level implementation
         * @param t timepoint
         * @param x Vector with the states
         * @param xdot Vector with the right hand side
         */
        void fxdot(realtype t, const_N_Vector x, N_Vector xdot);

        void fJ(realtype t, realtype cj, const Vector &x, const Vector &dx,
                const Vector &xdot, SUNMatrix J) override;

        /**
         * @brief Implementation of fJ at the N_Vector level
         *
         * This function provides an
         * interface to the model specific routines for the solver
         * implementation as well as the Vector level implementation
         * @param t timepoint
         * @param x Vector with the states
         * @param xdot Vector with the right hand side
         * @param J Matrix to which the Jacobian will be written
         **/
        void fJ(realtype t, const_N_Vector x, const_N_Vector xdot, SUNMatrix J);

        void fJSparse(realtype t, realtype cj, const Vector &x,
                      const Vector &dx, const Vector &xdot,
                      SUNMatrix J) override;

        /**
         * @brief Implementation of fJSparse at the N_Vector level, this function
         * provides an interface to the model specific routines for the solver
         * implementation as well as the Vector level implementation
         * @param t timepoint
         * @param x Vector with the states
         * @param J Matrix to which the Jacobian will be written
         */
        void fJSparse(realtype t, const_N_Vector x, SUNMatrix J);

        std::unique_ptr<Solver> getSolver() override;

    protected:

        /**
         * @brief Model specific implementation for fJSparse
         * @param JSparse Matrix to which the Jacobian will be written
         * @param t timepoint
         * @param x Vector with the states
         **/
        virtual void fJSparse(SUNMatrixContent_Sparse JSparse, realtype t, const realtype *x) = 0;

        /**
         * @brief Model specific implementation for fxdot
         * @param xdot residual function
         * @param t timepoint
         * @param x Vector with the states
         **/
        virtual void fxdot(realtype* xdot, realtype t, const realtype* x) = 0;

    };

}


#endif //SUNEIGEN_MODEL_ODE_H
