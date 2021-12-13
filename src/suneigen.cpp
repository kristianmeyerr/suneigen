#include "suneigen.h"

#include "misc.h"
#include "defines.h"
#include <cvodes/cvodes.h>           //return codes
#include <sundials/sundials_types.h> //realtype

#include <iostream>
#include <cstdarg>
#include <memory>

// ensure definitions are in sync
static_assert(suneigen::SUNEIGEN_SUCCESS == CV_SUCCESS,
              "SUNEIGEN_SUCCESS != CV_SUCCESS");
static_assert(suneigen::SUNEIGEN_DATA_RETURN == CV_TSTOP_RETURN,
              "SUNEIGEN_DATA_RETURN != CV_TSTOP_RETURN");
static_assert(suneigen::SUNEIGEN_ROOT_RETURN == CV_ROOT_RETURN,
              "SUNEIGEN_ROOT_RETURN != CV_ROOT_RETURN");
static_assert(suneigen::SUNEIGEN_ILL_INPUT == CV_ILL_INPUT,
              "SUNEIGEN_ILL_INPUT != CV_ILL_INPUT");
static_assert(suneigen::SUNEIGEN_NORMAL == CV_NORMAL,
              "SUNEIGEN_NORMAL != CV_NORMAL");
static_assert(suneigen::SUNEIGEN_ONE_STEP == CV_ONE_STEP,
              "SUNEIGEN_ONE_STEP != CV_ONE_STEP");
static_assert(suneigen::SUNEIGEN_SINGULAR_JACOBIAN == SUNLS_PACKAGE_FAIL_UNREC,
              "SUNEIGEN_SINGULAR_JACOBIAN != SUNLS_PACKAGE_FAIL_UNREC");
static_assert(std::is_same<suneigen::realtype, realtype>::value,
              "Definition of realtype does not match");

namespace suneigen {

    void printErrMsgIdAndTxt(std::string const& id, std::string const& message){
        std::cerr << "[Error] ";
        if (!id.empty()) {
            std::cerr << id << ": ";}
        std::cerr << message << std::endl;
    }

    void printWarnMsgIdAndTxt(std::string const& id, std::string const& message)
    {
        std::cerr << "[Warning] ";
        if (!id.empty()) {
            std::cerr << id << ": ";
        }
        std::cerr << message << std::endl;
    }

    std::unique_ptr<ReturnData>
    SunApplication::runSimulation(Solver& solver, Model& model, bool rethrow) {

        solver.startTimer();
        std::unique_ptr<ReturnData> rdata = std::make_unique<ReturnData>(solver, model);

        std::unique_ptr<ForwardProblem> fwd {};

        try {
            fwd = std::make_unique<ForwardProblem>(&model, &solver);

            fwd->workForwardProblem();

            rdata->status = SUNEIGEN_SUCCESS;

        } catch (suneigen::IntegrationFailure const& ex) {
            if(ex.error_code == SUNEIGEN_RHSFUNC_FAIL && solver.timeExceeded()) {
                rdata->status = SUNEIGEN_MAX_TIME_EXCEEDED;
                if(rethrow)
                    throw;
                warningF("SUNEIGEN:simulation",
                         "SUNEIGEN forward simulation failed at t = %f: "
                         "Maximum time exceeed.\n",
                         ex.time);
            } else {
                rdata->status = ex.error_code;
                if (rethrow)
                    throw;
                warningF("SUNEIGEN:simulation",
                         "SUNEIGEN forward simulation failed at t = %f:\n%s\n",
                         ex.time,
                         ex.what());

            }
        } catch (suneigen::SunException const& ex) {
            rdata->status = SUNEIGEN_ERROR;
            if (rethrow)
                throw;
            warningF("SUNEIGEN:simulation",
                     "SUNEIGEN simulation failed:\n%s\nError occurred in:\n%s",
                     ex.what(),
                     ex.getBacktrace());
        } catch (std::exception const& ex) {
            rdata->status = SUNEIGEN_ERROR;
            if (rethrow)
                throw;
            warningF("SUNEIGEN:simulation",
                     "SUNEIGEN simulation failed:\n%s\n",
                     ex.what());
        }

        rdata->processSimulationObjects(fwd.get(), model, solver);

        return rdata;
    }

    void SunApplication::warningF(const char* identifier, const char* format, ...) const
    {
        va_list argptr;
        va_start(argptr, format);
        auto str = printfToString(format, argptr);
        va_end(argptr);
        warning(identifier, str);
    }

    int SunApplication::checkFinite(gsl::span<const realtype> array, const char* fun)
    {

        for (unsigned int idx = 0; idx < array.size(); idx++) {
            if (std::isnan(array[idx])) {
                warningF("SUNEIGEN:NaN",
                         "SUNEIGEN encountered a NaN value at index %i/%i in %s!",
                         idx,
                         array.size()-1,
                         fun);
                return SUNEIGEN_RECOVERABLE_ERROR;
            }
            if (std::isinf(array[idx])) {
                warningF("SUNEIGEN:Inf",
                         "SUNEIGEN encountered an Inf value at index %i/%i in %s!",
                         idx,
                         array.size()-1,
                         fun);
                return SUNEIGEN_RECOVERABLE_ERROR;
            }
        }
        return SUNEIGEN_SUCCESS;
    }

}  // namespace suneigen
