
#ifndef SUNEIGEN_SOLVER_H
#define SUNEIGEN_SOLVER_H

#include <utility>
#include <functional>
#include <memory>
#include <vector>
#include <chrono>
#include <Eigen/Core>
#include <sundials/sundials_types.h>

#include "defines.h"

class Model;
class Solver;

class Solver {
public:

    Solver() = default;

    /**
     * @brief Solver copy constructor
     * @param other other.
     */
    Solver(const Solver &other);

    virtual ~Solver();

    /**
     * @brief Initializes the ami memory object and applies specified options
     * @param t0 initial timepoint
     * @param model pointer to the model instance
     * @param x0 initial states
     * @param dx0 initial derivative states
     * @param sx0 initial state sensitivities
     * @param sdx0 initial derivative state sensitivities
     */
     void setup(realtype t0, Model *model, const Eigen::VectorXd &x0,
             const Eigen::VectorXd &dx0, const Eigen::VectorXd &sx0,
             const Eigen::VectorXd &sdx0) const;

    /**
     * @brief Create specifies solver method and initializes solver memory for
     * the forward problem
     */
    virtual void allocateSolver() const = 0;


protected:

    /** pointer to solver memory block */
    mutable std::unique_ptr<void, std::function<void(void *)>> solver_memory_;

    /** specifies the linear multistep method.
     */
    LinearMultistepMethod lmm_ {LinearMultistepMethod::BDF};

private:

    /** absolute tolerances for integration */
    realtype atol_ {1e-16};


};

#endif //SUNEIGEN_SOLVER_H
