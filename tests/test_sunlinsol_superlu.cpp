#include <doctest/doctest.h>

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

#include <sunlinsol_superlu.h>

using namespace Eigen;

TEST_CASE("Test that we can solve a linear system for CSC_Mat"){

    int n = 4;

    // Fill matrix with uniform random data in [0,1/n]
    SUNMatrix B = SUNDenseMatrix(n, n);
    for (int k=0; k<5*n; k++) {
        int i = rand() % n;
        int j = rand() % n;
        realtype* matdata = SUNDenseMatrix_Column(B,j);
        matdata[i] = static_cast<realtype>(rand()) / static_cast<realtype>(RAND_MAX / n);
    }
    // Add identity to matrix
    int fails = SUNMatScaleAddI(1.0, B);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSol SUNMatScaleAddI failure\n");

    SUNMatrix A = SUNSparseFromDenseMatrix(B, 0.0, CSC_MAT);
    SUNMatDestroy(B);


    // Fill x vector with uniform random data in [0,1]
    N_Vector x = N_VNew_Serial(n);
    realtype* x_ptr = N_VGetArrayPointer(x);
    for (int i=0; i<n; i++)
        x_ptr[i] = static_cast<realtype>(rand()) / static_cast<realtype>(RAND_MAX);

    // create right-hand side vector for linear solve
    N_Vector b = N_VNew_Serial(n);
    fails += SUNMatMatvec(A, x, b);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSol SUNMatMatvec failure\n");


    // Create SuperLUMT linear solver
    SUNLinearSolver LS = SUNLinSol_SuperLU(x, A);

    SUNLinearSolver_Type mysuntype = SUNLinSolGetType(LS);
    REQUIRE(SUNLINEARSOLVER_DIRECT==mysuntype);

    sunindextype lastflag = SUNLinSolLastFlag(LS);
    REQUIRE(lastflag==0);

    fails += SUNLinSolInitialize(LS);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSolInitialize failure\n");
    fails += SUNLinSolSetup(LS, A);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSolSetup failure\n");
    fails += SUNLinSolSetZeroGuess(LS, true);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSolSetZeroGuess failure\n");

    N_Vector y = N_VNew_Serial(n);
    N_VScale(0.0, x, y);
    fails += SUNLinSolSolve(LS, A, y, b, -1);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSolSolve failure\n");

    Map<Matrix<double, Dynamic, 1> > x_map(N_VGetArrayPointer(x), n);
    Map<Matrix<double, Dynamic, 1> > y_map(N_VGetArrayPointer(y), n);

    REQUIRE_MESSAGE((x_map-y_map).sum() < 1e-14, "Fail: The solution error is to large ");

    sunindextype lastflag_after = SUNLinSolLastFlag(LS);
    REQUIRE(lastflag_after==SUNLS_SUCCESS);

    // Clean up
    SUNMatDestroy(A);
    N_VDestroy(x);
    N_VDestroy(b);
    SUNLinSolFree(LS);
    N_VDestroy(y);


}

TEST_CASE("Test that we can solve a linear system for CSR_MAT"){

    int n = 4;

    // Fill matrix with uniform random data in [0,1/n]
    SUNMatrix B = SUNDenseMatrix(n, n);
    for (int k=0; k<5*n; k++) {

        int i = rand() % n;
        int j = rand() % n;
        realtype* matdata = SUNDenseMatrix_Column(B,j);
        matdata[i] = static_cast<realtype>(rand()) / static_cast<realtype>(RAND_MAX / n);

    }
    // Add identity to matrix
    int fails = SUNMatScaleAddI(1.0, B);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSol SUNMatScaleAddI failure\n");

    SUNMatrix A = SUNSparseFromDenseMatrix(B, 0.0, CSR_MAT);
    SUNMatDestroy(B);

    // Fill x vector with uniform random data in [0,1]
    N_Vector x = N_VNew_Serial(n);
    realtype* x_ptr = N_VGetArrayPointer(x);
    for (int i=0; i<n; i++)
        x_ptr[i] = static_cast<realtype>(rand()) / static_cast<realtype>(RAND_MAX);

    // create right-hand side vector for linear solve
    N_Vector b = N_VNew_Serial(n);
    fails += SUNMatMatvec(A, x, b);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSol SUNMatMatvec failure\n");

    // Create SuperLUMT linear solver
    SUNLinearSolver LS = SUNLinSol_SuperLU(x, A);

    SUNLinearSolver_Type mysuntype = SUNLinSolGetType(LS);
    REQUIRE(SUNLINEARSOLVER_DIRECT==mysuntype);

    sunindextype lastflag = SUNLinSolLastFlag(LS);
    REQUIRE(lastflag==0);

    fails += SUNLinSolInitialize(LS);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSolInitialize failure\n");
    fails += SUNLinSolSetup(LS, A);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSolSetup failure\n");
    fails += SUNLinSolSetZeroGuess(LS, true);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSolSetZeroGuess failure\n");

    N_Vector y = N_VNew_Serial(n);
    N_VScale(0.0, x, y);
    fails += SUNLinSolSolve(LS, A, y, b, -1);
    REQUIRE_MESSAGE(fails==0, "FAIL: SUNLinSolSolve failure\n");

    Map<Matrix<double, Dynamic, 1> > x_map(N_VGetArrayPointer(x), n);
    Map<Matrix<double, Dynamic, 1> > y_map(N_VGetArrayPointer(y), n);

    REQUIRE_MESSAGE((x_map-y_map).sum() < 1e-14, "Fail: The solution error is to large ");

    sunindextype lastflag_after = SUNLinSolLastFlag(LS);
    REQUIRE(lastflag_after==SUNLS_SUCCESS);

    // Clean up
    SUNMatDestroy(A);
    N_VDestroy(b);
    N_VDestroy(x);
    N_VDestroy(y);
    SUNLinSolFree(LS);

}
