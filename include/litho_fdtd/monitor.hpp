#pragma once

#include <litho_fdtd/error.hpp>
#include <litho_fdtd/grid.hpp>
#include <litho_fdtd/spectrum.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace litho::fdtd {

// A monitor records the complex DFT of all six field components on an
// axis-aligned plane. `index` is the cell index along the plane normal.
struct monitor {
    std::string name;
    axis normal = axis::z;
    std::size_t index = 0;
    std::vector<double> frequencies; // empty means the source frequency
    std::optional<std::array<double, 2>> time_window; // seconds, inclusive
    std::optional<spectrum_options> spectrum;
};

struct monitor_result {
    std::string name;
    axis normal = axis::z;
    std::size_t index = 0;
    std::array<std::size_t, 2> shape{};
    std::array<double, 2> spacing{};
    std::array<double, 2> origin{};
    double normal_spacing = 0;
    std::vector<double> frequencies;
    std::vector<std::complex<double>> values; // [frequency][component][cell]
    std::optional<spectrum_result> spectrum;

    [[nodiscard]] std::size_t cell_count() const noexcept { return shape[0] * shape[1]; }
    [[nodiscard]] std::span<const std::complex<double>> component_values(std::size_t frequency_index,
                                                                        component value) const noexcept {
        const std::size_t cells = cell_count();
        const std::size_t channels = frequencies.size() * 6;
        const std::size_t channel = frequency_index * 6 + static_cast<std::size_t>(value);
        if (channel >= channels || cells == 0) return {};
        return {values.data() + channel * cells, cells};
    }
    [[nodiscard]] std::complex<double> at(std::size_t frequency_index, component value,
                                          std::size_t plane_index) const noexcept {
        const auto data = component_values(frequency_index, value);
        if (plane_index >= data.size()) return {};
        return data[plane_index];
    }
};

[[nodiscard]] inline result<void> validate(const monitor& value, const grid& domain) noexcept {
    const auto normal = static_cast<std::size_t>(value.normal);
    if (value.index >= domain.cells[normal]) return fail(errc::out_of_range);
    for (double frequency : value.frequencies)
        if (!std::isfinite(frequency) || !(frequency > 0)) return fail(errc::invalid_argument);
    if (value.time_window) {
        const auto& window = *value.time_window;
        if (!std::isfinite(window[0]) || !std::isfinite(window[1]) || !(window[0] < window[1]))
            return fail(errc::invalid_argument);
    }
    if (value.spectrum) {
        if (value.normal != axis::z) return fail(errc::unsupported);
        const auto& options = *value.spectrum;
        if (options.min_order_x > options.max_order_x || options.min_order_y > options.max_order_y)
            return fail(errc::invalid_argument);
        if (!value.frequencies.empty() && options.frequency_index >= value.frequencies.size())
            return fail(errc::invalid_argument);
        if (!std::isfinite(options.background_index.real()) || !(options.background_index.real() > 0)
            || !std::isfinite(options.background_index.imag()))
            return fail(errc::invalid_argument);
        if (!std::isfinite(options.period[0]) || !std::isfinite(options.period[1])
            || options.period[0] < 0 || options.period[1] < 0)
            return fail(errc::invalid_argument);
        if (options.propagation_sign != 1 && options.propagation_sign != -1)
            return fail(errc::invalid_argument);
        if (options.incident_power && (!std::isfinite(*options.incident_power) || *options.incident_power < 0))
            return fail(errc::invalid_argument);
        if (options.incident_wavevector) {
            const auto& wavevector = *options.incident_wavevector;
            if (!std::isfinite(wavevector[0]) || !std::isfinite(wavevector[1]))
                return fail(errc::invalid_argument);
        }
    }
    return {};
}

// Plane dimensions, spacing and the grid origin of the first sample.
[[nodiscard]] inline std::array<std::size_t, 2> plane_shape(const grid& domain, axis normal) noexcept {
    if (normal == axis::x) return {domain.cells[1], domain.cells[2]};
    if (normal == axis::y) return {domain.cells[0], domain.cells[2]};
    return {domain.cells[0], domain.cells[1]};
}

[[nodiscard]] inline std::array<double, 2> plane_spacing(const grid& domain, axis normal) noexcept {
    if (normal == axis::x) return {domain.spacing[1], domain.spacing[2]};
    if (normal == axis::y) return {domain.spacing[0], domain.spacing[2]};
    return {domain.spacing[0], domain.spacing[1]};
}

[[nodiscard]] inline std::array<double, 2> plane_origin(const grid& domain, axis normal) noexcept {
    if (normal == axis::x) return {domain.origin[1], domain.origin[2]};
    if (normal == axis::y) return {domain.origin[0], domain.origin[2]};
    return {domain.origin[0], domain.origin[1]};
}

} // namespace litho::fdtd
