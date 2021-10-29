//
// Created by Kristian Meyer on 2019-09-01.
//

#ifndef SUNEIGENMATRIX_NVECTOREIGEN_HPP
#define SUNEIGENMATRIX_NVECTOREIGEN_HPP

#include <sundials/sundials_nvector.h>
#include "Eigen/Dense"

namespace suneigen {

    struct _contentVectorXd {
        Eigen::VectorXd* vectorXdPointer;
        bool isOwnData;  // data ownership flag
    };

    // Forward declaration
    typedef struct _contentVectorXd* contentVectorXd;


    // --------------------------------------- //
    // Declare functions used by NVector_Eigen //
    // --------------------------------------- //

    /**
     * Allocates memory for a new NVector and attaches the operations and content structure with empty data
     *
     * @param [in] length Length of the NVector
     * @return Returns an empty NVector
     */
    N_Vector constructEmptyVectorXd(sunindextype length);

}

#endif //SUNEIGENMATRIX_NVECTOREIGEN_HPP
