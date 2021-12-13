#include "doctest/doctest.h"

#include "models/model_robertson.h"
#include "suneigen.h"

TEST_CASE("Test CVodes on Roberts problem using a dense linear solver."){

    // Create a model instance
    auto model = suneigen::generic_model::getModel();

    // Set desired output timepoints
    model->setTimepoints({4.0e1, 4.0e2, 4.0e8});

    // Set parameters
    std::vector<realtype> p{0.04, 1.0e4, 3.0e7};
    model->setParameters(p);

    // Set fixed parameters
    std::vector<realtype> k{1.0};
    model->setFixedParameters(k);

    // Create a solver instance
    auto solver = model->getSolver();

    // Optionally set integration tolerance
    solver->setAbsoluteTolerance(1e-12);
    solver->setRelativeTolerance(1e-8);
    solver->setLinearSolver(suneigen::LinearSolver::dense);

    // Create an application instance
    auto app = suneigen::SunApplication();

    // Run the simulation
    auto rdata = app.runSimulation(*solver,*model);

    // Expected values are taking from the SUNDIALS example cvsRoberts_dns.
    REQUIRE(doctest::Approx(rdata->x[0]).epsilon(0.0001) == 7.158017e-01);
    REQUIRE(doctest::Approx(rdata->x[1]).epsilon(0.000000001) == 9.1850e-06);
    REQUIRE(doctest::Approx(rdata->x[2]).epsilon(0.0001) == 2.8418e-01);
    REQUIRE(doctest::Approx(rdata->x[3]).epsilon(0.0001) == 4.5053e-01);
    REQUIRE(doctest::Approx(rdata->x[4]).epsilon(0.000000001) == 3.2232e-06);
    REQUIRE(doctest::Approx(rdata->x[5]).epsilon(0.0001) == 5.4946e-01);

    // The event output matches the time point of root in SUNDIALS example cvsRoberts_dns.
    REQUIRE(doctest::Approx(rdata->z[0]).epsilon(0.001) == 2.6391e-01);
    REQUIRE(doctest::Approx(rdata->z[1]).epsilon(0.001) == 2.0790e+07);

}

TEST_CASE("Test CVodes on Roberts problem using a sparse linear solver."){

    // Create a model instance
    auto model = suneigen::generic_model::getModel();

    // Set desired output timepoints
    model->setTimepoints({4.0e1, 4.0e2, 4.0e8});

    // Set parameters
    std::vector<realtype> p{0.04, 1.0e4, 3.0e7};
    model->setParameters(p);

    // Set fixed parameters
    std::vector<realtype> k{1.0};
    model->setFixedParameters(k);

    // Create a solver instance
    auto solver = model->getSolver();

    // Optionally set integration tolerance
    solver->setAbsoluteTolerance(1e-12);
    solver->setRelativeTolerance(1e-8);
    solver->setLinearSolver(suneigen::LinearSolver::SuperLU);

    // Create an application instance
    auto app = suneigen::SunApplication();

    // Run the simulation
    auto rdata = app.runSimulation(*solver,*model);

    // Expected values are taking from the SUNDIALS example called cvsRoberts_dns.
    REQUIRE(doctest::Approx(rdata->x[0]).epsilon(0.0001) == 7.158017e-01);
    REQUIRE(doctest::Approx(rdata->x[1]).epsilon(0.000000001) == 9.1850e-06);
    REQUIRE(doctest::Approx(rdata->x[2]).epsilon(0.0001) == 2.8418e-01);
    REQUIRE(doctest::Approx(rdata->x[3]).epsilon(0.0001) == 4.5053e-01);
    REQUIRE(doctest::Approx(rdata->x[4]).epsilon(0.000000001) == 3.2232e-06);
    REQUIRE(doctest::Approx(rdata->x[5]).epsilon(0.0001) == 5.4946e-01);

    // The event output matches the time point of root in SUNDIALS example cvsRoberts_dns.
    REQUIRE(doctest::Approx(rdata->z[0]).epsilon(0.001) == 2.6391e-01);
    REQUIRE(doctest::Approx(rdata->z[1]).epsilon(0.001) == 2.0790e+07);
}
