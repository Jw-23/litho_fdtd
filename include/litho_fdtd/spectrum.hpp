#pragma once

#include <litho_fdtd/error.hpp>
#include <litho_fdtd/grid.hpp>

#include <array>
#include <complex>
#include <cstddef>
#include <optional>
#include <vector>

namespace litho::fdtd {

struct monitor_result;

// Plane-wave decomposition of a z-normal monitor. Orders are indexed by
// (m, n) with kx = kx_incident + 2*pi*m/period_x, and the same convention in y.
// If `incident_wavevector` is unset the solver uses the source transverse
// wavevector, which is the correct choice for a periodic mask.
struct spectrum_options {
    int min_order_x = -4;
    int max_order_x = 4;
    int min_order_y = -4;
    int max_order_y = 4;
    std::size_t frequency_index = 0;
    std::array<double, 2> period{}; // 0 means the full grid extent
    std::complex<double> background_index{1.0, 0.0};
    std::optional<std::array<double, 2>> incident_wavevector; // rad/m
    int propagation_sign = 1; // +1 for +z transmission, -1 for -z reflection
    std::optional<double> incident_power; // W/m^2, enables efficiencies
};

struct diffraction_order {
    int m = 0;
    int n = 0;
    double kx = 0;
    double ky = 0;
    std::complex<double> kz{};
    std::array<std::complex<double>, 3> field{}; // Ex, Ey, Ez
    bool propagating = false;
    double power = 0; // time-averaged z Poynting flux, W/m^2
    std::optional<double> efficiency;
};

struct spectrum_result {
    double frequency = 0;
    std::vector<diffraction_order> orders;
    double total_power = 0;
    std::optional<double> incident_power;
    std::optional<double> total_efficiency;

    [[nodiscard]] const diffraction_order* find(int m, int n) const noexcept {
        for (const auto& order : orders)
            if (order.m == m && order.n == n) return &order;
        return nullptr;
    }
};

[[nodiscard]] result<spectrum_result> compute_spectrum(const monitor_result& monitor,
                                                      const spectrum_options& options) noexcept;

// The (0,0) order power is the natural incident reference for a plane-wave run.
[[nodiscard]] inline std::optional<double> zero_order_power(const spectrum_result& value) noexcept {
    if (const auto* order = value.find(0, 0)) return order->power;
    return std::nullopt;
}

// Divide every complex order by the dominant component of the reference (0,0)
// order. This yields a mask spectrum normalized to the incident wave while
// preserving the relative phase of all orders. Powers are left untouched.
[[nodiscard]] result<spectrum_result> normalize_spectrum(const spectrum_result& value,
                                                        const spectrum_result& reference) noexcept;

} // namespace litho::fdtd
