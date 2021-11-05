//
// Created by Kristian Meyer on 04/11/2021.
//

#ifndef SUNEIGEN_SOLVER_CVODES_H
#define SUNEIGEN_SOLVER_CVODES_H

#include "solver.h"

/**
 * @brief The CVodeSolver class is a wrapper around the SUNDIALS CVODES solver.
 */
class CVodeSolver : public Solver {
public:
    using Solver::Solver;

    ~CVodeSolver() override = default;

    void allocateSolver() const override;

};
#endif //SUNEIGEN_SOLVER_CVODES_H
