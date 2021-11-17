#include <iostream>

#include "model_robertson.h"
#include "suneigen.h"

int main() {

    // Create a model instance
    auto model = suneigen::generic_model::getModel();

    // Set desired output timepoints
    model->setTimepoints({4.0e1, 4.0e2});

    // Create a solver instance
    auto solver = model->getSolver();

    // Optionally set integration tolerance
    solver->setAbsoluteTolerance(1e-16);
    solver->setRelativeTolerance(1e-8);

    // Create an application instance
    auto app = suneigen::SunApplication();

    // Run the simulation
    auto rdata = app.runSimulation(*solver,*model);

    for (unsigned int i = 0; i < model->getTimepoints().size(); ++i) {
        std::cout << rdata->x[i] << " " << rdata->x[i + 1] << " " << rdata->x[i + 2] << std::endl;
    }

    // Statistics
    int numnonlinsolvconvfails_total = 0;
    for (auto& n : rdata->numnonlinsolvconvfails)
        numnonlinsolvconvfails_total += n;

    int numerrtestfails_total = 0;
    for (auto& n : rdata->numerrtestfails)
        numerrtestfails_total += n;

    int numrhsevals_total = 0;
    for (auto& n : rdata->numrhsevals)
        numrhsevals_total += n;

    int numsteps_total = 0;
    for (auto& n : rdata->numsteps)
        numsteps_total += n;

    std::cout << "CPU time (s): " << rdata->cpu_time << std::endl;
    std::cout << "Number of RHS evaluations: " << numrhsevals_total << std::endl;
    std::cout << "Number of steps used: " << numsteps_total << std::endl;
    std::cout << "Number of nonlinear solver convergence failures: " << numnonlinsolvconvfails_total << std::endl;
    std::cout << "Number of error test failures: " << numerrtestfails_total << std::endl;

    /*

    const int n = 3;
    SUNMatrix A_dense = SUNDenseMatrix(n, n);
    realtype* A_dense_ptr = SUNDenseMatrix_Data(A_dense);
    A_dense_ptr[0] = 1;
    A_dense_ptr[2] = 7;
    A_dense_ptr[4] = 3;
    A_dense_ptr[8] = 5;
    SUNMatrix A = SUNSparseFromDenseMatrix(A_dense, 0.0, CSR_MAT);
    SUNSparseMatrix_Print(A, stdout);

    N_Vector b = N_VNew_Serial(n);
    realtype* b_ptr = N_VGetArrayPointer(b);
    b_ptr[0] = 1;
    b_ptr[1] = 3;
    b_ptr[2] = 2;

    N_Vector x = N_VNew_Serial(n);
    realtype* x_ptr = N_VGetArrayPointer(x);

    // Map N_Vectors
    Map<Matrix<double, n, 1> > x_map(N_VGetArrayPointer(x), n);
    Map<Matrix<double, n, 1> > b_map(N_VGetArrayPointer(b), n);

    // Map SUNMat
    sunindextype rows = SUNSparseMatrix_Rows(A);
    sunindextype cols = SUNSparseMatrix_Columns(A);
    sunindextype nnz = SUNSparseMatrix_NNZ(A);
    sunindextype* outerIndexPtr = SUNSparseMatrix_IndexPointers(A);
    sunindextype* innerIndicesPtr = SUNSparseMatrix_IndexValues(A);
    realtype* values = SUNSparseMatrix_Data(A);

    if(SM_SPARSETYPE_S(A)==CSC_MAT){
        Map<const SparseMatrix<double, ColMajor, sunindextype> > A_map(rows,cols,nnz,outerIndexPtr,innerIndicesPtr,values);
        SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.analyzePattern(A_map);
        solver.factorize(A_map);
        x_map = solver.solve(b_map);
    } else {
        Map<const SparseMatrix<double, RowMajor, sunindextype> > A_map(rows,cols,nnz,outerIndexPtr,innerIndicesPtr,values);
        SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.analyzePattern(A_map);
        solver.factorize(A_map);
        x_map = solver.solve(b_map);
    }

    // The solution is now stored in the N_Vector x
    for (int i = 0; i < n; ++i) {
        cout << x_ptr[i] << endl;
    }

    cout << typeid(sunindextype).name() << endl;

    */

    /*
    // fill b
    // solve Ax = b
    SolverClassName<SparseMatrix<double> > solver;
    solver.compute(A);
    if(solver.info()!=Success) {
        // decomposition failed
        return;
    }
    x = solver.solve(b);
    if(solver.info()!=Success) {
        // solving failed
        return;
    }
// solve for another right hand side:
    x1 = solver.solve(b1);

    // SUNDIALS solve
    const long int cols = 3;
    N_Vector b = N_VNew_Serial(cols);
    realtype* b_ptr = N_VGetArrayPointer(b);
    b_ptr[0] = 3.;
    b_ptr[1] = 3.;
    b_ptr[2] = 4.;

    SUNMatrix A = SUNDenseMatrix(cols, cols);
    realtype* A_ptr = SUNDenseMatrix_Data(A);
    A_ptr[0] = 1;
    A_ptr[1] = 4;
    A_ptr[2] = 7;
    A_ptr[3] = 2;
    A_ptr[4] = 5;
    A_ptr[5] = 8;
    A_ptr[6] = 3;
    A_ptr[7] = 6;
    A_ptr[8] = 10;
    SUNDenseMatrix_Print(A,stdout);

    N_Vector x = N_VNew_Serial(cols);
    N_VConst(0.0, x);

    SUNMatrix AA = SUNDenseMatrix(cols, cols);
    SUNMatCopy(A, AA);

    int fails = 0;
    SUNLinearSolver LS = SUNLinSol_Dense(x, A);
    LS->ops->solve = partial_piv_lu;
    LS->ops->setup = nullptr;

    fails = SUNLinSolInitialize(LS);
    fails = SUNLinSolSetup(LS, A);
    fails = SUNLinSolSetZeroGuess(LS, true);
    fails = LS->ops->solve(LS, A, x, b, 1e-14);

    N_VPrint(x);
    SUNDenseMatrix_Print(A,stdout);

    // Eigen map SUNDIALS vec/matrix solve
    Map<Matrix<double, cols, 1> > b_map(N_VGetArrayPointer(b), cols);
    Map<Matrix<double, cols, cols> > A_map(SUNDenseMatrix_Data(AA), cols, cols);
    PartialPivLU<Ref<MatrixXd> > lu(A_map);
    VectorXd x_map = lu.solve(b_map);

    cout << "The solution is:\n" << x_map << endl;
    cout << "A is:\n" << A_map << endl;
    */
    return 0;
}
