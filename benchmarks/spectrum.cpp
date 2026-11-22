#include <litho_fdtd/fdtd.hpp>
#include <chrono>
#include <complex>
#include <cstdio>
#include <numbers>

int main() {
    using namespace litho::fdtd;
    monitor_result plane;
    plane.shape = {120, 120};
    plane.spacing = {5e-9, 5e-9};
    plane.origin = {13e-9, -17e-9};
    plane.frequencies = {speed_of_light / 193e-9};
    plane.values.resize(6 * plane.cell_count());
    for (std::size_t c = 0; c < 6; ++c)
        for (std::size_t j = 0; j < 120; ++j)
            for (std::size_t i = 0; i < 120; ++i)
                plane.values[c * plane.cell_count() + j * 120 + i] =
                    std::polar(1.0 / (c + 1), 2 * std::numbers::pi * (2.0 * i - j) / 120);
    spectrum_options options;
    options.min_order_x = options.min_order_y = -6;
    options.max_order_x = options.max_order_y = 6;
    constexpr int repeats = 5;
    const auto start = std::chrono::steady_clock::now();
    double checksum = 0;
    for (int i = 0; i < repeats; ++i) {
        auto result = compute_spectrum(plane, options);
        if (!result) return 1;
        checksum += std::abs(result->find(2, -1)->field[0]);
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("spectrum: 120 x 120, 169 orders, %.6f s/extraction, checksum %.12g\n", seconds / repeats, checksum);
}
