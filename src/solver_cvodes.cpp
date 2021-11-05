#include "solver_cvodes.h"

#include <cvodes/cvodes.h>

void CVodeSolver::allocateSolver() const {
    if (!solver_memory_)
        solver_memory_ = std::unique_ptr<void, std::function<void(void *)>>(
                CVodeCreate(static_cast<int>(lmm_)),
                [](void *ptr) { CVodeFree(&ptr); });
}
