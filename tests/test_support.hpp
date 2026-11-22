#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>

namespace test {

inline int failures = 0;

inline bool close(double a, double b, double tolerance) noexcept {
    if (!std::isfinite(a) || !std::isfinite(b)) return false;
    const double scale = std::max({1.0e-12, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= tolerance * scale;
}

inline bool close(std::complex<double> a, std::complex<double> b, double tolerance) noexcept {
    return close(a.real(), b.real(), tolerance) && close(a.imag(), b.imag(), tolerance);
}

inline void report(const char* file, int line, const char* expression) {
    std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", file, line, expression);
    ++failures;
}

} // namespace test

#define CHECK(expression)                                                                          \
    do {                                                                                           \
        if (!(expression)) test::report(__FILE__, __LINE__, #expression);                          \
    } while (false)

#define CHECK_NEAR(actual, expected, tolerance)                                                    \
    do {                                                                                           \
        if (!test::close((actual), (expected), (tolerance))) {                                     \
            std::fprintf(stderr, "%s:%d: CHECK_NEAR failed: %s = %.12g, expected %.12g\n",         \
                         __FILE__, __LINE__, #actual, static_cast<double>(actual),                 \
                         static_cast<double>(expected));                                           \
            ++test::failures;                                                                      \
        }                                                                                          \
    } while (false)

#define CHECK_COMPLEX_NEAR(actual, expected, tolerance)                                            \
    do {                                                                                           \
        if (!test::close((actual), (expected), (tolerance))) {                                     \
            std::fprintf(stderr, "%s:%d: CHECK_COMPLEX_NEAR failed: %s = (%.12g, %.12g), "          \
                         "expected (%.12g, %.12g)\n",                                              \
                         __FILE__, __LINE__, #actual, (actual).real(), (actual).imag(),            \
                         (expected).real(), (expected).imag());                                    \
            ++test::failures;                                                                      \
        }                                                                                          \
    } while (false)
