
#ifndef SUNEIGEN_SUNLINSOL_SUPERLU_H
#define SUNEIGEN_SUNLINSOL_SUPERLU_H

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_sparse.h>

#include "Eigen/SparseLU"

// SuperLU Implementation of SUNLinearSolver
struct _SUNLinearSolverContent_SuperLU {
    int last_flag;
    int  first_factorize;
    sunindextype n;
    Eigen::SparseLU<Eigen::SparseMatrix<double> >* solver;
};
typedef struct _SUNLinearSolverContent_SuperLU* SUNLinearSolverContent_SuperLU;

/**
 * The function SUNLinSol_SuperLU creates and allocates memory for a SuperLU
 * SUNLinearSolver object.
 *
 * @param [in] y Input vector
 * @param [in] A Input Matrix
 * @return Returns the SUNLinearSolver object
 */
SUNDIALS_EXPORT SUNLinearSolver SUNLinSol_SuperLU(N_Vector y, SUNMatrix A);

/**
 * The required function SUNLinSolGetType returns the type identifier for the
 * linear solver LS. It is used to determine the solver type
 * (direct, iterative, or matrix-iterative) from the abstract SUNLinearSolver interface.
 *
 * @param [in] S Input linear solver
 * @return Returns the Solver type
 */
SUNDIALS_EXPORT SUNLinearSolver_Type SUNLinSolGetType_SuperLU(SUNLinearSolver S);

/**
 * The required function SUNLinSolInitialize_SuperLU performs linear solver initialization
 * (assuming that all solver-specific options have been set).
 *
 * @param [in] S Input linear solver
 * @return Returns zero for a successful call, and a negative value for a failure,
 *      ideally returning one of the generic error codes
 */
SUNDIALS_EXPORT int SUNLinSolInitialize_SuperLU(SUNLinearSolver S);

/**
 * The first time that the “setup” routine is called,
 * it performs the symbolic factorization, followed by an initial numerical factorization.
 * On subsequent calls to the “setup” routine, it skips the symbolic factorization,
 * and only refactors the input matrix.
 *
 * @param [in] S Input linear solver
 * @param [in] A Input matrix
 * @return Returns success/failure
 */
SUNDIALS_EXPORT int SUNLinSolSetup_SuperLU(SUNLinearSolver S, SUNMatrix A);

/**
 * The required function SUNLinSolSolve solves a linear system @f$ Ax = b @f$.
 *
 * @param [in] S Input sunlinsol object.
 * @param [in] A A is sunmatrix object.
 * @param [out] x x is a nvector object containing the initial guess for the solution of the linear system,
 *      and the solution to the linear system upon return.
 * @param [in] b b is an nvector object containing the linear system right-hand side.
 * @param [in] tol Tol is the desired linear solver tolerance (not used in direct solver)
 * @return This should return zero for a successful call, a positive value
 *      for a recoverable failure and a negative value for an unrecoverable failure,
 *      ideally returning one of the generic error codes.
 */
SUNDIALS_EXPORT int SUNLinSolSolve_SuperLU(SUNLinearSolver S, SUNMatrix A,
                                             N_Vector x, N_Vector b, realtype tol);

SUNDIALS_EXPORT sunindextype SUNLinSolLastFlag_SuperLU(SUNLinearSolver S);

SUNDIALS_EXPORT int SUNLinSolSpace_SuperLU(SUNLinearSolver S,
                                             long int *lenrwLS,
                                             long int *leniwLS);

/**
 * Frees memory allocated by the linear solver
 * @param [in] S Input Linear solver
 * @return Returns success/failure
 */
SUNDIALS_EXPORT int SUNLinSolFree_SuperLU(SUNLinearSolver S);

// Macro-like functions
int& FIRSTFACTORIZE(SUNLinearSolver S);
int& LASTFLAG(SUNLinearSolver S);
Eigen::SparseLU<Eigen::SparseMatrix<double> >& SOLVER(SUNLinearSolver S);

#endif //SUNEIGEN_SUNLINSOL_SUPERLU_H
