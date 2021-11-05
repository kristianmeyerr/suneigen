#include "suneigen_exception.h"
#include "misc.h"


#include <cstdarg>
#include <cstdio>
#include <cstring>

SunException::SunException()
{
    backtraceString();
}

SunException::SunException(const char *fmt, ...)
        : SunException()
{
    va_list ap;
    va_start(ap, fmt);
    storeMessage(fmt, ap);
    va_end(ap);
}

const char *SunException::what() const noexcept {
    return msg_.data();
}

const char *SunException::getBacktrace() const {
    return trace_.data();
}
void SunException::storeBacktrace() {
    snprintf(trace_.data(), trace_.size(), "%s",
             backtraceString().c_str());
}

// https://stackoverflow.com/questions/20167124/vsprintf-and-vsnprintf-wformat-nonliteral-warning-on-clang-5-0
__attribute__ ((__format__ (__printf__, 2, 0)))
void SunException::storeMessage(const char *fmt, va_list argptr)
{
    vsnprintf(msg_.data(), msg_.size(), fmt, argptr);
}

CvodeException::CvodeException(const int error_code, const char *function) :
        SunException("Cvode routine %s failed with error code %i",function,error_code){}

__attribute__ ((__format__ (__printf__, 2, 0)))
SetupFailure::SetupFailure(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    storeMessage(fmt, ap);
    va_end(ap);
}

void SunException::outOfLineVirtualFix() {}
void CvodeException::outOfLineVirtualFixCvode() {}
void SetupFailure::outOfLineVirtualFixSetup() {}
