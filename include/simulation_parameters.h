#ifndef SUNEIGEN_SIMULATION_PARAMETERS_H
#define SUNEIGEN_SIMULATION_PARAMETERS_H

#include "defines.h"

#include <iostream>

namespace suneigen {

    /**
     * @brief Container for various simulation parameters.
     */
    class SimulationParameters {
    public:
        SimulationParameters() = default;

        /**
         * @brief Constructor
         * @param fixedParameters_ Model constants
         * @param parameters_ Model parameters
         */
        SimulationParameters(std::vector<realtype> fixedParameters_,
                             std::vector<realtype> parameters_)
                : fixedParameters(std::move(fixedParameters_)),
                  parameters(std::move(parameters_))
        {}

        /**
         * @brief Model constants
         *
         * Vector of size Model::nk() or empty
         */
        std::vector<realtype> fixedParameters;

        /**
         * @brief Model parameters
         *
         * Vector of size Model::np() or empty
         */
        std::vector<realtype> parameters;

        /**
         * @brief Initial state
         *
         * Vector of size Model::nx() or empty
         */
        std::vector<realtype> x0;

        /**
         * @brief Initial state sensitivities
         *
         * Dimensions:
         * Model::nx() * Model::nplist(), Model::nx() * ExpData::plist.size(), if
         * ExpData::plist is not empty, or empty
         */
        std::vector<realtype> sx0;

        /** starting time */
        realtype tstart_{0.0};

        /**
         * @brief Timepoints for which model state/outputs/... are requested
         *
         * Vector of timepoints.
         */
        std::vector<realtype> ts_;
    };

    bool operator==(const SimulationParameters &a, const SimulationParameters &b);

}

#endif //SUNEIGEN_SIMULATION_PARAMETERS_H
