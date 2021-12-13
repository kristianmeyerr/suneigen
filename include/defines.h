#ifndef SUNEIGEN_DEFINES_H
#define SUNEIGEN_DEFINES_H

#include <functional>
#include <string>
#include <cmath>

namespace suneigen {

    /* Return codes */
    constexpr int SUNEIGEN_RECOVERABLE_ERROR = 1;
    constexpr int SUNEIGEN_UNRECOVERABLE_ERROR = -10;
    constexpr int SUNEIGEN_TOO_MUCH_WORK = -1;
    constexpr int SUNEIGEN_TOO_MUCH_ACC = -2;
    constexpr int SUNEIGEN_ERR_FAILURE = -3;
    constexpr int SUNEIGEN_CONV_FAILURE = -4;
    constexpr int SUNEIGEN_RHSFUNC_FAIL = -8;
    constexpr int SUNEIGEN_ILL_INPUT = -22;
    constexpr int SUNEIGEN_ERROR = -99;
    constexpr int SUNEIGEN_NOT_IMPLEMENTED = -999;
    constexpr int SUNEIGEN_MAX_TIME_EXCEEDED = -1000;
    constexpr int SUNEIGEN_SINGULAR_JACOBIAN = -809;
    constexpr int SUNEIGEN_SUCCESS = 0;
    constexpr int SUNEIGEN_DATA_RETURN = 1;
    constexpr int SUNEIGEN_ROOT_RETURN = 2;

    constexpr int SUNEIGEN_NORMAL = 1;
    constexpr int SUNEIGEN_ONE_STEP = 2;

    /** defines variable type for simulation variables
     * (determines numerical accuracy) */
    using realtype = double;

    /** methods for sensitivity computation */
    enum class SensitivityMethod {
        none,
        forward,
        adjoint
    };

    /** CVODES/IDAS forward sensitivity computation method */
    enum class InternalSensitivityMethod {
        simultaneous = 1,
        staggered = 2,
        staggered1 = 3
    };

    /** CVODES/IDAS state interpolation for adjoint sensitivity analysis */
    enum class InterpolationType {
        hermite = 1,
        polynomial = 2
    };

    /** CVODES/IDAS Nonlinear Iteration method */
    enum class NonlinearSolverIteration {
        fixedpoint = 1,
        newton = 2
    };

    /** linear solvers for CVODES/IDAS */
    enum class LinearSolver {
        dense = 1,
        SuperLU = 2,
    };

    /** orders of sensitivity analysis */
    enum class SensitivityOrder {
        none,
        first,
        second
    };

    /** CVODES/IDAS linear multistep method */
    enum class LinearMultistepMethod {
        adams = 1,
        BDF = 2
    };

    /** Damping factor flag for the Newton method */
    enum class NewtonDampingFactorMode {
        off = 0,
        on = 1
    };

    /**
     * Type for function to process warnings or error messages.
     */
    using outputFunctionType = std::function<void(std::string const &identifier,
                                                  std::string const &message)>;

}  // namespace suneigen

#endif //SUNEIGEN_DEFINES_H
