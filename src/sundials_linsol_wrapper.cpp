#include <sundials_linsol_wrapper.h>

#include <exception.h>

#include <new> // bad_alloc
#include <utility>

namespace suneigen {

    SUNLinSolWrapper::SUNLinSolWrapper(SUNLinearSolver linsol) : solver_(linsol) {}

    SUNLinSolWrapper::~SUNLinSolWrapper() {
        if (solver_)
            SUNLinSolFree(solver_);
    }

    SUNLinSolWrapper::SUNLinSolWrapper(SUNLinSolWrapper &&other) noexcept {
        std::swap(solver_, other.solver_);
    }

    SUNLinearSolver SUNLinSolWrapper::get() const { return solver_; }

    SUNLinearSolver_Type SUNLinSolWrapper::getType() const {
        return SUNLinSolGetType(solver_);
    }

    int SUNLinSolWrapper::initialize() {
        auto res = SUNLinSolInitialize(solver_);
        if (res != SUNLS_SUCCESS)
            throw SunException("Solver initialization failed with code %d", res);
        return res;
    }

    void SUNLinSolWrapper::setup(SUNMatrix A) const {
        auto res = SUNLinSolSetup(solver_, A);
        if (res != SUNLS_SUCCESS)
            throw SunException("Solver setup failed with code %d", res);
    }

    void SUNLinSolWrapper::setup(const SUNMatrixWrapper& A) const { return setup(A.get()); }

    int SUNLinSolWrapper::Solve(SUNMatrix A, N_Vector x, N_Vector b, realtype tol) const {
        return SUNLinSolSolve(solver_, A, x, b, tol);
    }

    long SUNLinSolWrapper::getLastFlag() const {
        return static_cast<long>(SUNLinSolLastFlag(solver_));
    }

    int SUNLinSolWrapper::space(long *lenrwLS, long *leniwLS) const {
        return SUNLinSolSpace(solver_, lenrwLS, leniwLS);
    }

    SUNMatrix SUNLinSolWrapper::getMatrix() const { return nullptr; }

    SUNNonLinSolWrapper::SUNNonLinSolWrapper(SUNNonlinearSolver sol)
            : solver(sol) {}

    SUNNonLinSolWrapper::~SUNNonLinSolWrapper() {
        if (solver)
            SUNNonlinSolFree(solver);
    }

    SUNNonLinSolWrapper::SUNNonLinSolWrapper(SUNNonLinSolWrapper &&other) noexcept {
        std::swap(solver, other.solver);
    }

    SUNNonLinSolWrapper &SUNNonLinSolWrapper::
    operator=(SUNNonLinSolWrapper &&other) noexcept {
        std::swap(solver, other.solver);
        return *this;
    }

    SUNNonlinearSolver SUNNonLinSolWrapper::get() const { return solver; }

    SUNNonlinearSolver_Type SUNNonLinSolWrapper::getType() const {
        return SUNNonlinSolGetType(solver);
    }

    int SUNNonLinSolWrapper::setup(N_Vector y, void *mem) {
        auto res = SUNNonlinSolSetup(solver, y, mem);
        if (res != SUN_NLS_SUCCESS)
            throw SunException("Nonlinear solver setup failed with code %d", res);
        return res;
    }

    int SUNNonLinSolWrapper::Solve(N_Vector y0, N_Vector y, N_Vector w,
                                   realtype tol, bool callLSetup, void *mem) {
        return SUNNonlinSolSolve(solver, y0, y, w, tol, callLSetup, mem);
    }

    int SUNNonLinSolWrapper::setSysFn(SUNNonlinSolSysFn SysFn) {
        return SUNNonlinSolSetSysFn(solver, SysFn);
    }

    int SUNNonLinSolWrapper::setLSetupFn(SUNNonlinSolLSetupFn SetupFn) {
        return SUNNonlinSolSetLSetupFn(solver, SetupFn);
    }

    int SUNNonLinSolWrapper::setLSolveFn(SUNNonlinSolLSolveFn SolveFn) {
        return SUNNonlinSolSetLSolveFn(solver, SolveFn);
    }

    int SUNNonLinSolWrapper::setConvTestFn(SUNNonlinSolConvTestFn CTestFn,
                                           void* ctest_data) {
        return SUNNonlinSolSetConvTestFn(solver, CTestFn, ctest_data);
    }

    int SUNNonLinSolWrapper::setMaxIters(int maxiters) {
        return SUNNonlinSolSetMaxIters(solver, maxiters);
    }

    long SUNNonLinSolWrapper::getNumIters() const {
        long int niters = -1;
        auto res = SUNNonlinSolGetNumIters(solver, &niters);
        if (res != SUN_NLS_SUCCESS) {
            throw SunException("SUNNonlinSolGetNumIters failed with code %d", res);
        }
        return niters;
    }

    int SUNNonLinSolWrapper::getCurIter() const {
        int iter = -1;
        auto res = SUNNonlinSolGetCurIter(solver, &iter);
        if (res != SUN_NLS_SUCCESS) {
            throw SunException("SUNNonlinSolGetCurIter failed with code %d", res);
        }
        return iter;
    }

    long SUNNonLinSolWrapper::getNumConvFails() const {
        long int nconvfails = -1;
        auto res = SUNNonlinSolGetNumConvFails(solver, &nconvfails);
        if (res != SUN_NLS_SUCCESS) {
            throw SunException("SUNNonlinSolGetNumConvFails failed with code %d",
                               res);
        }
        return nconvfails;
    }

    void SUNNonLinSolWrapper::initialize() {
        int status = SUNNonlinSolInitialize(solver);
        if (status != SUN_NLS_SUCCESS)
            throw SunException(
                    "Nonlinear solver initialization failed with code %d", status);
    }

    SUNLinSolDense::SUNLinSolDense(const Vector &x) :
    A_(SUNMatrixWrapper(static_cast<sunindextype>(x.getLength()),
            static_cast<sunindextype>(x.getLength())))
    {
        solver_ = SUNLinSol_Dense(const_cast<N_Vector>(x.getNVector()), A_.get());
        if (!solver_)
            throw SunException("Failed to create solver.");
    }

    SUNMatrix SUNLinSolDense::getMatrix() const { return A_.get(); }

    SUNLinSolSuperLU::SUNLinSolSuperLU(N_Vector x_, SUNMatrix A_)
            : SUNLinSolWrapper(SUNLinSol_SuperLU(x_, A_))
    {
        if (!solver_)
            throw SunException("Failed to create solver.");
    }

    SUNLinSolSuperLU::SUNLinSolSuperLU(Vector &x, int nnz, int sparsetype)
            : A(SUNMatrixWrapper(static_cast<sunindextype>(x.getLength()),
                    static_cast<sunindextype>(x.getLength()), nnz, sparsetype))
    {
        solver_ = SUNLinSol_SuperLU(x.getNVector(), A.get());
        if (!solver_)
            throw SunException("Failed to create solver.");
    }

    SUNMatrix SUNLinSolSuperLU::getMatrix() const
    {
        return A.get();
    }


    SUNNonLinSolNewton::SUNNonLinSolNewton(N_Vector x)
            : SUNNonLinSolWrapper(SUNNonlinSol_Newton(x)) {
    }

    SUNNonLinSolNewton::SUNNonLinSolNewton(int count, N_Vector x)
            : SUNNonLinSolWrapper(SUNNonlinSol_NewtonSens(count, x)) {
        if (!solver)
            throw(SunException("SUNNonlinSol_NewtonSens failed"));
    }

    int SUNNonLinSolNewton::getSysFn(SUNNonlinSolSysFn *SysFn) const {
        return SUNNonlinSolGetSysFn_Newton(solver, SysFn);
    }

    SUNNonLinSolFixedPoint::SUNNonLinSolFixedPoint(const_N_Vector x, int m)
            : SUNNonLinSolWrapper(SUNNonlinSol_FixedPoint(const_cast<N_Vector>(x), m)) {
    }

    SUNNonLinSolFixedPoint::SUNNonLinSolFixedPoint(int count, const_N_Vector x, int m)
            : SUNNonLinSolWrapper(
            SUNNonlinSol_FixedPointSens(count, const_cast<N_Vector>(x), m)) {
    }

    int SUNNonLinSolFixedPoint::getSysFn(SUNNonlinSolSysFn *SysFn) const {
        return SUNNonlinSolGetSysFn_FixedPoint(solver, SysFn);
    }

} // namespace suneigen
