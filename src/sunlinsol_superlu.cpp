//
// Created by Kristian Meyer on 03/11/2021.
//

#include "sunlinsol_superlu.h"
#include <sundials/sundials_math.h>
#include <iostream>

SUNLinearSolver SUNLinSol_SuperLU(N_Vector y, SUNMatrix A){

    // Require sparse and square matrix
    if (SUNMatGetID(A) != SUNMATRIX_SPARSE) return(nullptr);
    if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A)) return(nullptr);

    if ( (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
         (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
         (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS) )
        return(nullptr);

    sunindextype MatrixRows = SUNSparseMatrix_Rows(A);
    if (MatrixRows != N_VGetLength(y)) return(nullptr);

    // Create an empty linear solver
    SUNLinearSolver S = SUNLinSolNewEmpty();
    if (S == nullptr) return(nullptr);

    // Attach operations
    S->ops->gettype    = SUNLinSolGetType_SuperLU;
    S->ops->getid      = nullptr;
    S->ops->initialize = SUNLinSolInitialize_SuperLU;
    S->ops->setup      = SUNLinSolSetup_SuperLU;
    S->ops->solve      = SUNLinSolSolve_SuperLU;
    S->ops->lastflag   = SUNLinSolLastFlag_SuperLU;
    S->ops->space      = nullptr;
    S->ops->free       = SUNLinSolFree_SuperLU;

    // Create content
    SUNLinearSolverContent_SuperLU content = nullptr;
    content = reinterpret_cast<SUNLinearSolverContent_SuperLU>(malloc(sizeof *content));
    if (content == nullptr) { SUNLinSolFree(S); return(nullptr); }

    // Attach content
    S->content = content;

    // Fill content
    content->n = MatrixRows;
    content->last_flag = 0;

    // Use new since it calls the Eigen constructor
    try {
        content->solver = new Eigen::SparseLU<Eigen::SparseMatrix<double> >;
    } catch(const std::bad_alloc& e){
        SUNLinSolFree(S);
        std::cout << "Allocation failed: " << e.what() << std::endl;
        return(nullptr);
    }

    return(S);
}

SUNLinearSolver_Type SUNLinSolGetType_SuperLU(SUNLinearSolver S){
    (void)S;
    return(SUNLINEARSOLVER_DIRECT);
}

int SUNLinSolInitialize_SuperLU(SUNLinearSolver S)
{
    /* force a first factorization */
    FIRSTFACTORIZE(S) = 1;

    LASTFLAG(S) = SUNLS_SUCCESS;
    return(LASTFLAG(S));
}

int SUNLinSolSetup_SuperLU(SUNLinearSolver S, SUNMatrix A){

    // Get an Eigen map to A
    sunindextype rows = SUNSparseMatrix_Rows(A);
    sunindextype cols = SUNSparseMatrix_Columns(A);
    sunindextype nnz = SUNSparseMatrix_NNZ(A);
    sunindextype* outerIndexPtr = SUNSparseMatrix_IndexPointers(A);
    sunindextype* innerIndicesPtr = SUNSparseMatrix_IndexValues(A);
    realtype* values = SUNSparseMatrix_Data(A);

    if(SM_SPARSETYPE_S(A)==CSC_MAT){

        // Column major map
        Eigen::Map<Eigen::SparseMatrix<double, Eigen::ColMajor, sunindextype> > A_map(rows,cols,nnz,outerIndexPtr,innerIndicesPtr,values);

        // On first decomposition, set up reusable pieces
        if (FIRSTFACTORIZE(S)) {
            // Perform symbolic factorization
            SOLVER(S).analyzePattern(A_map);
            FIRSTFACTORIZE(S) = 0;
        }

        // Compute the LU factorization of A.
        SOLVER(S).factorize(A_map);
        if (SOLVER(S).info() != Eigen::Success) {
            LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
            return LASTFLAG(S);
        }

    } else {

        // Row major map
        Eigen::Map<Eigen::SparseMatrix<double, Eigen::RowMajor, sunindextype> > A_map(rows,cols,nnz,outerIndexPtr,innerIndicesPtr,values);

        // On first decomposition, set up reusable pieces
        if (FIRSTFACTORIZE(S)) {
            // Perform symbolic factorization
            SOLVER(S).analyzePattern(A_map);
            FIRSTFACTORIZE(S) = 0;
        }

        // Compute the LU factorization of A.
        SOLVER(S).factorize(A_map);
        if (SOLVER(S).info() != Eigen::Success) {
            LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
            return LASTFLAG(S);
        }
    }

    LASTFLAG(S) = SUNLS_SUCCESS;
    return(LASTFLAG(S));
}

int SUNLinSolSolve_SuperLU(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                             N_Vector b, realtype tol){

    // tol is un-used by direct solver
    (void) tol;

    sunindextype rows = SUNSparseMatrix_Rows(A);
    sunindextype cols = SUNSparseMatrix_Columns(A);
    sunindextype nnz = SUNSparseMatrix_NNZ(A);
    sunindextype* outerIndexPtr = SUNSparseMatrix_IndexPointers(A);
    sunindextype* innerIndicesPtr = SUNSparseMatrix_IndexValues(A);
    realtype* values = SUNSparseMatrix_Data(A);

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> > x_map(N_VGetArrayPointer(x), rows);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 1> > b_map(N_VGetArrayPointer(b), rows);

    if(SM_SPARSETYPE_S(A)==CSC_MAT){
        Eigen::Map<const Eigen::SparseMatrix<double, Eigen::ColMajor, sunindextype> > A_map(rows,cols,nnz,outerIndexPtr,innerIndicesPtr,values);
        x_map = SOLVER(S).solve(b_map);
        if (SOLVER(S).info() != Eigen::Success) {
            LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
            return(LASTFLAG(S));
        }
    } else {
        Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, sunindextype> > A_map(rows,cols,nnz,outerIndexPtr,innerIndicesPtr,values);
        x_map = SOLVER(S).solve(b_map);
        if (SOLVER(S).info() != Eigen::Success) {
            LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
            return(LASTFLAG(S));
        }
    }

    LASTFLAG(S) = SUNLS_SUCCESS;
    return(LASTFLAG(S));
}

sunindextype SUNLinSolLastFlag_SuperLU(SUNLinearSolver S)
{
    /* return the stored 'last_flag' value */
    if (S == nullptr) return(-1);
    return(LASTFLAG(S));
}

int SUNLinSolFree_SuperLU(SUNLinearSolver S)
{
    // return with success if already freed
    if (S == nullptr) return(SUNLS_SUCCESS);

    // delete items from the contents structure (if it exists)
    if (S->content) {

        if(reinterpret_cast<SUNLinearSolverContent_SuperLU>(S->content)->solver){
            free(reinterpret_cast<SUNLinearSolverContent_SuperLU>(S->content)->solver);
            reinterpret_cast<SUNLinearSolverContent_SuperLU>(S->content)->solver = nullptr;
        }

        free(S->content);
        S->content = nullptr;
    }

    // delete generic structures
    if (S->ops) {
        free(S->ops);
        S->ops = nullptr;
    }
    free(S); S = nullptr;
    return(SUNLS_SUCCESS);
}

// Functions that replace macros in sundials style.
int& FIRSTFACTORIZE(SUNLinearSolver S)
{
  return reinterpret_cast<SUNLinearSolverContent_SuperLU>(S->content)->first_factorize;
}

int& LASTFLAG(SUNLinearSolver S)
{
  return reinterpret_cast<SUNLinearSolverContent_SuperLU>(S->content)->last_flag;
}

Eigen::SparseLU<Eigen::SparseMatrix<double> >& SOLVER(SUNLinearSolver S)
{
  return *reinterpret_cast<SUNLinearSolverContent_SuperLU>(S->content)->solver;
}
