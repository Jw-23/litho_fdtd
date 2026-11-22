#include <litho_fdtd/spectrum.hpp>

#include <litho_fdtd/constants.hpp>
#include <litho_fdtd/monitor.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>

namespace litho::fdtd {

namespace {

using complex = std::complex<double>;

} // namespace

result<spectrum_result> compute_spectrum(const monitor_result& monitor, const spectrum_options& options) noexcept {
    return detail::guard([&]() -> result<spectrum_result> {
        if (monitor.normal != axis::z) return fail(errc::unsupported);
        if (monitor.frequencies.empty()) return fail(errc::invalid_argument);
        if (options.frequency_index >= monitor.frequencies.size()) return fail(errc::out_of_range);
        if (options.min_order_x > options.max_order_x || options.min_order_y > options.max_order_y)
            return fail(errc::invalid_argument);
        if (options.incident_power && !(*options.incident_power > 0)) return fail(errc::invalid_argument);
        if (!std::isfinite(options.period[0]) || !std::isfinite(options.period[1]) ||
            options.period[0] < 0 || options.period[1] < 0)
            return fail(errc::invalid_argument);
        if (options.incident_wavevector) {
            const auto& wavevector = *options.incident_wavevector;
            if (!std::isfinite(wavevector[0]) || !std::isfinite(wavevector[1]))
                return fail(errc::invalid_argument);
        }

        const std::size_t nx = monitor.shape[0];
        const std::size_t ny = monitor.shape[1];
        if (nx == 0 || ny == 0) return fail(errc::invalid_argument);
        const double dx = monitor.spacing[0];
        const double dy = monitor.spacing[1];
        const double period_x = options.period[0] > 0 ? options.period[0] : static_cast<double>(nx) * dx;
        const double period_y = options.period[1] > 0 ? options.period[1] : static_cast<double>(ny) * dy;
        if (!(period_x > 0) || !(period_y > 0)) return fail(errc::invalid_argument);

        const double frequency = monitor.frequencies[options.frequency_index];
        const double omega = 2.0 * std::numbers::pi * frequency;
        const double k0 = omega / speed_of_light;
        const complex k = options.background_index * k0;
        const double kx_incident = options.incident_wavevector ? (*options.incident_wavevector)[0] : 0.0;
        const double ky_incident = options.incident_wavevector ? (*options.incident_wavevector)[1] : 0.0;

        const auto ex = monitor.component_values(options.frequency_index, component::ex);
        const auto ey = monitor.component_values(options.frequency_index, component::ey);
        const auto ez = monitor.component_values(options.frequency_index, component::ez);
        const auto hx = monitor.component_values(options.frequency_index, component::hx);
        const auto hy = monitor.component_values(options.frequency_index, component::hy);
        if (ex.size() != nx * ny || ey.size() != nx * ny || ez.size() != nx * ny ||
            hx.size() != nx * ny || hy.size() != nx * ny)
            return fail(errc::invalid_configuration);

        spectrum_result result;
        result.frequency = frequency;
        result.incident_power = options.incident_power;
        const std::size_t order_count_x = static_cast<std::size_t>(options.max_order_x - options.min_order_x + 1);
        const std::size_t order_count_y = static_cast<std::size_t>(options.max_order_y - options.min_order_y + 1);
        result.orders.reserve(order_count_x * order_count_y);

        std::vector<complex> phase_x(nx);
        for (int n = options.min_order_y; n <= options.max_order_y; ++n) {
            const double ky = ky_incident + 2.0 * std::numbers::pi * static_cast<double>(n) / period_y;
            for (int m = options.min_order_x; m <= options.max_order_x; ++m) {
                const double kx = kx_incident + 2.0 * std::numbers::pi * static_cast<double>(m) / period_x;

                complex sum_ex{};
                complex sum_ey{};
                complex sum_ez{};
                complex sum_hx{};
                complex sum_hy{};
                // Separate x/y phases and apply Yee offsets once per order.
                // No transcendental functions remain in the cell loop.
                for (std::size_t i = 0; i < nx; ++i)
                    phase_x[i] = std::exp(complex{0.0, -kx * (monitor.origin[0] + i * dx)});
                for (std::size_t j = 0; j < ny; ++j) {
                    const complex phase_y = std::exp(complex{0.0, -ky * (monitor.origin[1] + j * dy)});
                    for (std::size_t i = 0; i < nx; ++i) {
                        const std::size_t plane_index = i + nx * j;
                        const complex phase = phase_x[i] * phase_y;
                        sum_ex += ex[plane_index] * phase;
                        sum_ey += ey[plane_index] * phase;
                        sum_ez += ez[plane_index] * phase;
                        sum_hx += hx[plane_index] * phase;
                        sum_hy += hy[plane_index] * phase;
                    }
                }
                const complex offset_x = std::exp(complex{0.0, -0.5 * kx * dx});
                const complex offset_y = std::exp(complex{0.0, -0.5 * ky * dy});
                sum_ex *= offset_x;
                sum_ey *= offset_y;
                sum_ez *= offset_x * offset_y;
                sum_hx *= offset_y;
                sum_hy *= offset_x;
                const double norm = 1.0 / static_cast<double>(nx * ny);
                const std::array<complex, 3> field{sum_ex * norm, sum_ey * norm, sum_ez * norm};

                const complex kz_squared = k * k - complex{kx * kx + ky * ky, 0.0};
                complex kz = std::sqrt(kz_squared);
                if (kz.imag() < 0) kz = -kz;
                if (options.propagation_sign < 0) kz = -kz;

                const bool propagating = std::abs(kz.real()) > 1e-12 * std::abs(k) &&
                                          std::abs(kz.real()) >= std::abs(kz.imag());

                double power = 0;
                if (propagating) {
                    // H is sampled half a cell above the E plane. Remove that
                    // spatial phase before forming the Poynting flux.
                    const complex stagger = std::exp(complex{0.0, -1.0} * kz * (0.5 * monitor.normal_spacing));
                    const complex hx_t = sum_hx * norm * stagger;
                    const complex hy_t = sum_hy * norm * stagger;
                    const complex sz = 0.5 * (field[0] * std::conj(hy_t) - field[1] * std::conj(hx_t));
                    power = static_cast<double>(options.propagation_sign) * sz.real();
                }

                diffraction_order order;
                order.m = m;
                order.n = n;
                order.kx = kx;
                order.ky = ky;
                order.kz = kz;
                order.field = field;
                order.propagating = propagating;
                order.power = power;
                if (propagating && options.incident_power && *options.incident_power > 0)
                    order.efficiency = power / *options.incident_power;
                result.total_power += power;
                result.orders.push_back(order);
            }
        }

        if (options.incident_power && *options.incident_power > 0)
            result.total_efficiency = result.total_power / *options.incident_power;
        return result;
    });
}

result<spectrum_result> normalize_spectrum(const spectrum_result& value,
                                          const spectrum_result& reference) noexcept {
    return detail::guard([&]() -> result<spectrum_result> {
        const diffraction_order* incident = reference.find(0, 0);
        if (!incident) return fail(errc::invalid_argument);
        std::size_t dominant = 0;
        double magnitude = std::abs(incident->field[0]);
        for (std::size_t c = 1; c < 3; ++c) {
            const double candidate = std::abs(incident->field[c]);
            if (candidate > magnitude) {
                magnitude = candidate;
                dominant = c;
            }
        }
        const complex scale = incident->field[dominant];
        if (!(magnitude > 0)) return fail(errc::numerical_failure);
        spectrum_result out = value;
        for (auto& order : out.orders)
            for (auto& component_value : order.field) component_value /= scale;
        return out;
    });
}

} // namespace litho::fdtd
