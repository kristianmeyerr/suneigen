#ifndef SUNEIGEN_MISC_H
#define SUNEIGEN_MISC_H

#include "vector.h"
#include "exception.h"

#include <gsl/gsl-lite.hpp>
#include <string>

namespace suneigen {

    /**
     * @brief creates a slice from existing data
     *
     * @param data to be sliced
     * @param index slice index
     * @param size slice size
     * @return span of the slice
     */
    template <class T>
    gsl::span<T> slice(std::vector<T> &data, size_t index, size_t size) {
        if ((index + 1) * size > data.size())
            throw std::out_of_range("requested slice is out of data range");
        if (size > 0)
            return gsl::make_span(&data.at(index*size), size);

        return gsl::make_span(static_cast<T*>(nullptr), 0);
    }

    /**
     * @brief creates a constant slice from existing constant data
     *
     * @param data to be sliced
     * @param index slice index
     * @param size slice size
     * @return span of the slice
     */
    template <class T>
    gsl::span<const T> slice(const std::vector<T> &data,
                             size_t index, size_t size) {
        if ((index + 1) * size > data.size())
            throw std::out_of_range("requested slice is out of data range");
        if (size > 0)
            return gsl::make_span(&data.at(index*size), size);

        return gsl::make_span(static_cast<T*>(nullptr), 0);
    }

    /**
     * @brief local helper to check whether the provided buffer has the expected
     * size
     * @param buffer buffer to which values are to be written
     * @param expected_size expected size of the buffer
     */
    template <class T>
    void checkBufferSize(gsl::span<T> buffer,
                         typename gsl::span<T>::index_type expected_size) {
        if (buffer.size() != expected_size)
            throw SunException("Incorrect buffer size! Was %u, expected %u.",
                               buffer.size(), expected_size);
    }

    /**
     * @brief local helper function to write computed slice to provided buffer (span)
     * @param slice computed value
     * @param buffer buffer to which values are to be written
     */
    template <class T>
    void writeSlice(const gsl::span<const T> slice, gsl::span<T> buffer) {
        checkBufferSize(buffer, slice.size());
        std::copy(slice.begin(), slice.end(), buffer.data());
    }

    /**
     * @brief local helper function to write computed slice to provided buffer (vector)
     * @param s computed value
     * @param b buffer to which values are to be written
     */
    template <class T>
    void writeSlice(const std::vector<T> &s, std::vector<T> &b) {
        writeSlice(gsl::make_span(s.data(), s.size()),
                   gsl::make_span(b.data(), b.size()));
    }

    /**
     * @brief local helper function to write computed slice to provided buffer (vector/span)
     * @param s computed value
     * @param b buffer to which values are to be written
     */
    template <class T>
    void writeSlice(const std::vector<T> &s, gsl::span<T> b) {
        writeSlice(gsl::make_span(s.data(), s.size()), b);
    }

    /**
     * @brief local helper function to write computed slice to provided buffer (Vector/span)
     * @param s computed value
     * @param b buffer to which values are to be written
     */
    void writeSlice(const Vector &s, gsl::span<realtype> b);

    /**
     * @brief Returns the current backtrace as std::string
     * @param skip Number of frames to skip
     * @return Backtrace
     */
    std::string backtraceString(int skip = 2);

    /**
     * @brief Format printf-style arguments to std::string
     * @param fmt Format string
     * @param ap Argument list pointer
     * @return Formatted String
     */
    std::string printfToString(const char *fmt, va_list ap);

    /**
     * @brief Generic implementation for a context manager, explicitly deletes copy
     * and move operators for derived classes
     */
    class ContextManager{
    public:
        ContextManager() = default;
        ContextManager(ContextManager &other) = delete;
        ContextManager(ContextManager &&other) = delete;
    };

}

#endif //SUNEIGEN_MISC_H
