#include "doctest/doctest.h"
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include "nvector_eigen.hpp"

int check_ans(realtype ans, N_Vector X, sunindextype local_length);
int check_ans(realtype ans, N_Vector X, sunindextype local_length){

    realtype* Xdata = N_VGetArrayPointer(X);

    // check vector data
    int failure = 0;
    for (int i = 0; i < local_length; i++) {
        failure += SUNRCompare(Xdata[i], ans);
    }

    return (failure > 0) ? (1) : (0);
}

TEST_CASE("Test that we can create NVectors"){
    N_Vector W = N_VNewEmpty_Eigen(100);
    REQUIRE(W != nullptr);

    N_Vector X = N_VNew_Eigen(100);
    REQUIRE(X != nullptr);
    REQUIRE(N_VGetArrayPointer(X) != nullptr);

    N_VDestroy(W);
    N_VDestroy(X);
}

TEST_CASE("Check vector ID")
{
    N_Vector X = N_VNew_Eigen(100);
    REQUIRE(N_VGetVectorID(X) == 15);
    N_VDestroy(X);
}

TEST_CASE("Check vector length"){
    N_Vector X = N_VNew_Eigen(100);

    // check if the required operations are implemented
    REQUIRE(X->ops->nvconst != nullptr);
    REQUIRE(X->ops->nvdotprod != nullptr);

    sunindextype n1 = N_VGetLength(X);
    N_VConst(RCONST(1.0), X);

    // use N_VConst and N_VDotProd to compute length
    N_VConst(RCONST(1.0), X);
    auto n2 = static_cast<sunindextype>(N_VDotProd(X, X));

    REQUIRE(n1==n2);

    N_VDestroy(X);
}

TEST_CASE("Test clone functions"){

    int len = 100;
    N_Vector W = N_VNew_Eigen(len);

    // Clone empty vector
    N_Vector X = N_VCloneEmpty(W);
    REQUIRE(X != nullptr);

    // Clone vector with data
    N_Vector Z = N_VClone(W);
    REQUIRE(Z != nullptr);
    REQUIRE(N_VGetArrayPointer(Z) != nullptr);

    // check if the required operations are implemented
    REQUIRE(Z->ops->nvconst != nullptr);

    N_VConst(1.0,Z);
    REQUIRE(check_ans(1.0, Z, len) == 0);

    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Z);

}




