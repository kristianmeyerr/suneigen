#include "doctest/doctest.h"

#include "model_robertson.h"
#include "suneigen.h"

TEST_CASE("Test CVodes on Roberts problem using a dense Jacobian."){

    // Create a model instance
    auto model = suneigen::generic_model::getModel();

    // Set desired output timepoints
    model->setTimepoints({4.0e1, 4.0e2});

    // Create a solver instance
    auto solver = model->getSolver();

    // Optionally set integration tolerance
    solver->setAbsoluteTolerance(1e-12);
    solver->setRelativeTolerance(1e-8);

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
}
