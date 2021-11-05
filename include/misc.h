#ifndef SUNEIGEN_MISC_H
#define SUNEIGEN_MISC_H


#include <string>

/**
 * @brief Returns the current backtrace as std::string
 * @param skip Number of frames to skip
 * @return Backtrace
 */
std::string backtraceString(int skip = 2);

#endif //SUNEIGEN_MISC_H
