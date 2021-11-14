#ifndef SUNEIGEN_SUNEIGEN_H
#define SUNEIGEN_SUNEIGEN_H

#include "defines.h"
#include "exception.h"
#include "model.h"
#include "solver.h"
#include "return_data.h"


namespace suneigen {

    /*!
     * @brief Prints a specified error message associated with the specified
     * identifier
     *
     * @param id error identifier
     * @param message error message
     */
    void printErrMsgIdAndTxt(std::string const &id, std::string const &message);

    /*!
     * @brief Prints a specified warning message associated with the specified
     * identifier
     *
     * @param id warning identifier
     * @param message warning message
     */
    void printWarnMsgIdAndTxt(std::string const &id, std::string const &message);

    /**
     * @brief Main class for making calls to SUNEIGEN.
     *
     * This class is used to provide separate SUNEIGEN contexts, for example, for use
     * in multi-threaded applications where different threads want to use SUNEIGEN with
     * different settings, such custom logging functions.
     *
     * NOTE: For this moment, the context object needs to be set manually to any
     * Model and Solver object. If not set, they will use the default output
     * channel.
     */
    class SunApplication {
    public:
        SunApplication() = default;

        /**
         * @brief Core integration routine. Initializes the solver and runs the
         * forward and backward problem.
         *
         * @param solver Solver instance
         * @param model model specification object
         * @param rethrow rethrow integration exceptions?
         * @return rdata pointer to return data object
         */
        std::unique_ptr<ReturnData> runSimulation(Solver &solver,Model &model, bool rethrow = false);

        /** Function to process warnings */
        outputFunctionType warning = printWarnMsgIdAndTxt;

        /** Function to process errors */
        outputFunctionType error = printErrMsgIdAndTxt;

        /**
         * @brief printf interface to warning()
         * @param identifier warning identifier
         * @param format string with warning message printf-style format
         * @param ... arguments to be formatted
         */
        void warningF(const char *identifier, const char *format, ...) const;

        /**
         * @brief Checks the values in an array for NaNs and Infs
         *
         * @param array array
         * @param fun name of calling function
         * @return AMICI_RECOVERABLE_ERROR if a NaN/Inf value was found,
         * AMICI_SUCCESS otherwise
         */
        int checkFinite(gsl::span<const realtype> array, const char *fun);
    };

}

#endif //SUNEIGEN_SUNEIGEN_H
