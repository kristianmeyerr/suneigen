#ifndef SUNEIGEN_EXCEPTION_H
#define SUNEIGEN_EXCEPTION_H

#include <exception>
#include <array>

#include "defines.h"

namespace suneigen {

    /**
     * @brief Suneigen exception class
     *
     * Has a printf style interface to allow easy generation of error messages
     */
    class SunException : public std::exception {
    public:
        /**
         * @brief Constructor storing backtrace
         */
        SunException();

        /**
         * @brief Constructor with printf style interface
         * @param fmt error message with printf format
         * @param ... printf formatting variables
         */
        explicit SunException(char const* fmt, ...);

        /**
         * @brief Override of default error message function
         * @return msg error message
         */
        [[nodiscard]] const char* what() const noexcept override;

        /**
         * @brief Returns the stored backtrace
         * @return trace backtrace
         */
        [[nodiscard]] const char *getBacktrace() const;

        /**
         * @brief Stores the current backtrace
         */
        void storeBacktrace();

    protected:
        /**
         * @brief Store the provided message
         * @param fmt error message with printf format
         * @param argptr pointer to variadic argument list
         */
        void storeMessage(const char *fmt, va_list argptr);

    private:

        std::array<char, 500> msg_{};
        std::array<char, 500> trace_{};

    };

    /**
     * @brief cvode exception handler class
     */
    class CvodeException : public SunException  {
    public:
        /**
         * @brief Constructor
         * @param error_code error code returned by cvode function
         * @param function cvode function name
         */
        CvodeException(int error_code, const char *function);

    };

    /**
     * @brief Integration failure exception for the forward problem
     *
     * This exception should be thrown when an integration failure occurred
     * for this exception we can assume that we can recover from the exception
     * and return a solution struct to the user
     */
    class IntegrationFailure : public SunException  {
    public:
        /**
         * @brief Constructor
         * @param code error code returned by cvode/ida
         * @param t time of integration failure
         */
        IntegrationFailure(int code, realtype t);

        /** error code returned by cvodes/idas */
        int error_code;

        /** time of integration failure */
        realtype time;
    };

    /**
     * @brief Setup failure exception
     *
     * This exception should be thrown when the solver setup failed
     * for this exception we can assume that we cannot recover from the exception
     * and an error will be thrown
     */
    class SetupFailure : public SunException {
    public:
        /**
         * @brief Constructor with printf style interface
         * @param fmt error message with printf format
         * @param ... printf formatting variables
         */
        explicit SetupFailure(char const* fmt, ...);


    };

}

#endif //SUNEIGEN_EXCEPTION_H
