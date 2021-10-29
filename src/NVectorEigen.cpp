//
// Created by Kristian Meyer on 2019-09-01.
//

#include <sundials/sundials_math.h>
#include "NVectorEigen.hpp"
#include <iostream>
#include "Eigen/Dense"

namespace suneigen {

    static int constructEmptyVectorXd(int length)
    {
        return length;
    }
}
