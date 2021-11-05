#ifndef SUNEIGEN_DEFINES_H
#define SUNEIGEN_DEFINES_H

/** linear solvers for CVODES/IDAS */
enum class LinearSolver {
    dense     = 1,
    band      = 2,
    SuperLU   = 3,
};

/** CVODES/IDAS linear multistep method */
enum class LinearMultistepMethod {
    adams = 1,
    BDF = 2
};

#endif //SUNEIGEN_DEFINES_H
