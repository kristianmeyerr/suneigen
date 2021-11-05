//
// Created by Kristian Meyer on 04/11/2021.
//

#include "solver.h"
#include "suneigen_exception.h"

Solver::~Solver() = default;

Solver::Solver(const Solver &other) : atol_(other.atol_) {}

void Solver::setup(realtype t0, Model *model, const Eigen::VectorXd &x0,
           const Eigen::VectorXd &dx0, const Eigen::VectorXd &sx0,
           const Eigen::VectorXd &sdx0) const {
    (void)t0;
    (void)model;
    (void)x0;
    (void)dx0;
    (void)sx0;
    (void)sdx0;

    /* Create solver memory object if necessary */
    allocateSolver();
    if (!solver_memory_)
        throw SunException("Failed to allocated solver memory!");
}
