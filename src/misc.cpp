// https://gist.github.com/fmela/591333

#include "misc.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>

#if defined(_WIN32)
#define PLATFORM_WINDOWS // Windows
#elif defined(_WIN64)
#define PLATFORM_WINDOWS // Windows
#elif defined(__CYGWIN__) && !defined(_WIN32)
#define PLATFORM_WINDOWS // Windows (Cygwin POSIX under Microsoft Window)
#else
#include <execinfo.h>
#include <dlfcn.h>    // for dladdr
#include <cxxabi.h>   // for __cxa_demangle
#endif


// This function produces a stack backtrace with demangled function & method names.
std::string backtraceString(int skip)
{
    std::ostringstream trace_buf;

#ifdef PLATFORM_WINDOWS
    trace_buf << "stacktrace not available on windows platforms\n";
#else
    void *callstack[128];
    const int nMaxFrames = sizeof(callstack) / sizeof(callstack[0]);
    char buf[1024];
    int nFrames = backtrace(callstack, nMaxFrames);
    char **symbols = backtrace_symbols(callstack, nFrames);

    // Skip the first to omit SunEigenException and storeBacktrace
    for (int i = skip; i < nFrames; i++) {
        // call
        Dl_info info;
        if (dladdr(callstack[i], &info) && info.dli_sname) {
            char *demangled = nullptr;
            int status = -1;
            if (info.dli_sname[0] == '_')
                demangled = abi::__cxa_demangle(info.dli_sname, nullptr, nullptr,
                                                &status);
            snprintf(buf, sizeof(buf), "%-3d %*p %s + %zd\n", i - 2,
                     int(2 + sizeof(void *) * 2), callstack[i],
                     status == 0 ? demangled
                                 : info.dli_sname == nullptr ? symbols[i]
                                                             : info.dli_sname,
                     static_cast<ssize_t>(static_cast<char *>(callstack[i]) -
                                          static_cast<char *>(info.dli_saddr)));

            free(demangled);
        } else {
            snprintf(buf, sizeof(buf), "%-3d %*p %s\n", i - 2,
                     int(2 + sizeof(void *) * 2), callstack[i],
                     symbols[i]);
        }
        trace_buf << buf;
    }
    free(symbols);

    if (nFrames == nMaxFrames)
        trace_buf << "[truncated]\n";
#endif
    return trace_buf.str();
}
