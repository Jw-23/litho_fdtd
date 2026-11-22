#include <litho_fdtd/solver.hpp>

#include "detail/cpml.hpp"
#include "detail/parallel_for.hpp"

#include <litho_fdtd/spectrum.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <type_traits>
#include <utility>
#include <vector>

namespace litho::fdtd {

namespace {

using complex = std::complex<double>;

template<class Scalar>
struct field_state {
    std::vector<Scalar> ex, ey, ez;
    std::vector<Scalar> hx, hy, hz;
};

template<class Scalar>
void allocate(field_state<Scalar>& fields, std::size_t count) {
    fields.ex.assign(count, Scalar{});
    fields.ey.assign(count, Scalar{});
    fields.ez.assign(count, Scalar{});
    fields.hx.assign(count, Scalar{});
    fields.hy.assign(count, Scalar{});
    fields.hz.assign(count, Scalar{});
}

template<class Scalar>
struct dispersion_state {
    std::vector<Scalar> jx, jy, jz;
};

struct material_coefficients {
    std::array<std::vector<double>, 3> ca; // Ex, Ey, Ez
    std::array<std::vector<double>, 3> cb;
    std::array<std::vector<double>, 3> h_cb; // Hx, Hy, Hz
    std::array<std::vector<double>, 3> drude_a;
    std::array<std::vector<double>, 3> drude_b;
    bool has_drude = false;
};

struct source_info {
    double frequency = 0;
    double omega = 0;
    double kx = 0, ky = 0, kz = 0, k_magnitude = 0;
    double z0 = 0;
    double sigma_z = 1;
    double pulse_sigma = 1;
    double phase_offset = 0;
    std::array<double, 3> e0{};
    std::array<double, 3> h0{};
    std::array<double, 3> h1{}; // z-hat cross e0 / (omega mu)
    std::array<complex, 3> phase{complex{1}, complex{1}, complex{1}};
    std::array<complex, 3> phase_inv{complex{1}, complex{1}, complex{1}};
    bool complex_mode = false;
};

template<class Scalar>
Scalar apply_phase(Scalar value, complex phase) noexcept {
    if constexpr (std::is_same_v<Scalar, double>)
        return value * phase.real();
    else
        return value * phase;
}

template<class Scalar>
Scalar diff_backward(const Scalar* data, std::size_t index, std::size_t stride, std::size_t coordinate,
                     std::size_t count, bool periodic, complex phase_inverse, double inverse_spacing) noexcept {
    if (coordinate > 0) return (data[index] - data[index - stride]) * inverse_spacing;
    if (periodic) {
        const Scalar wrapped = apply_phase(data[index + (count - 1) * stride], phase_inverse);
        return (data[index] - wrapped) * inverse_spacing;
    }
    return (data[index + stride] - data[index]) * inverse_spacing;
}

template<class Scalar>
Scalar diff_forward(const Scalar* data, std::size_t index, std::size_t stride, std::size_t coordinate,
                    std::size_t count, bool periodic, complex phase, double inverse_spacing) noexcept {
    if (coordinate + 1 < count) return (data[index + stride] - data[index]) * inverse_spacing;
    if (periodic) {
        const Scalar wrapped = apply_phase(data[index - (count - 1) * stride], phase);
        return (wrapped - data[index]) * inverse_spacing;
    }
    return (data[index] - data[index - stride]) * inverse_spacing;
}

template<class Scalar>
result<void> update_ex(const grid& domain, field_state<Scalar>& f, const detail::cpml_set& cpml,
               detail::cpml_memory<Scalar>& memory, const material_coefficients& coefficients,
               dispersion_state<Scalar>* dispersion, complex phase_y_inverse, complex phase_z_inverse,
               std::size_t workers) {
    const std::size_t nx = domain.cells[0], ny = domain.cells[1], nz = domain.cells[2];
    const std::size_t sy = nx, sz = nx * ny;
    const double inverse_dy = 1.0 / domain.spacing[1], inverse_dz = 1.0 / domain.spacing[2];
    const auto& axis_y = cpml.axes[1];
    const auto& axis_z = cpml.axes[2];
    Scalar* ex = f.ex.data();
    const Scalar* hy = f.hy.data();
    const Scalar* hz = f.hz.data();
    Scalar* psi_y = memory.ex_y.empty() ? nullptr : memory.ex_y.data();
    Scalar* psi_z = memory.ex_z.empty() ? nullptr : memory.ex_z.data();
    const auto* ca = coefficients.ca[0].data();
    const auto* cb = coefficients.cb[0].data();
    const auto* drude_a = coefficients.has_drude ? coefficients.drude_a[0].data() : nullptr;
    const auto* drude_b = coefficients.has_drude ? coefficients.drude_b[0].data() : nullptr;

    return detail::parallel_for(nz * ny, workers, [&](std::size_t row_begin, std::size_t row_end) {
        for (std::size_t row = row_begin; row < row_end; ++row) {
            const std::size_t k = row / ny;
            const std::size_t j = row % ny;
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t index = i + nx * (j + ny * k);
                const Scalar dhy_dz =
                    diff_backward(hy, index, sz, k, nz, axis_z.periodic, phase_z_inverse, inverse_dz);
                const Scalar dhz_dy =
                    diff_backward(hz, index, sy, j, ny, axis_y.periodic, phase_y_inverse, inverse_dy);
                Scalar stretched_y = dhz_dy;
                Scalar stretched_z = dhy_dz;
                if (!axis_y.periodic) {
                    const Scalar previous = psi_y[index];
                    stretched_y = dhz_dy * axis_y.inv_kappa_int[j] + previous;
                    psi_y[index] = axis_y.b_e_int[j] * previous + axis_y.c_e_int[j] * dhz_dy;
                }
                if (!axis_z.periodic) {
                    const Scalar previous = psi_z[index];
                    stretched_z = dhy_dz * axis_z.inv_kappa_int[k] + previous;
                    psi_z[index] = axis_z.b_e_int[k] * previous + axis_z.c_e_int[k] * dhy_dz;
                }
                const Scalar curl = stretched_y - stretched_z;
                if (dispersion) {
                    const Scalar previous = dispersion->jx[index];
                    const Scalar next = drude_a[index] * previous + drude_b[index] * ex[index];
                    ex[index] = ca[index] * ex[index] + cb[index] * (curl - 0.5 * (next + previous));
                    dispersion->jx[index] = next;
                } else {
                    ex[index] = ca[index] * ex[index] + cb[index] * curl;
                }
            }
        }
    });
}

template<class Scalar>
result<void> update_ey(const grid& domain, field_state<Scalar>& f, const detail::cpml_set& cpml,
               detail::cpml_memory<Scalar>& memory, const material_coefficients& coefficients,
               dispersion_state<Scalar>* dispersion, complex phase_x_inverse, complex phase_z_inverse,
               std::size_t workers) {
    const std::size_t nx = domain.cells[0], ny = domain.cells[1], nz = domain.cells[2];
    const std::size_t sx = 1, sz = nx * ny;
    const double inverse_dx = 1.0 / domain.spacing[0], inverse_dz = 1.0 / domain.spacing[2];
    const auto& axis_x = cpml.axes[0];
    const auto& axis_z = cpml.axes[2];
    Scalar* ey = f.ey.data();
    const Scalar* hx = f.hx.data();
    const Scalar* hz = f.hz.data();
    Scalar* psi_x = memory.ey_x.empty() ? nullptr : memory.ey_x.data();
    Scalar* psi_z = memory.ey_z.empty() ? nullptr : memory.ey_z.data();
    const auto* ca = coefficients.ca[1].data();
    const auto* cb = coefficients.cb[1].data();
    const auto* drude_a = coefficients.has_drude ? coefficients.drude_a[1].data() : nullptr;
    const auto* drude_b = coefficients.has_drude ? coefficients.drude_b[1].data() : nullptr;

    return detail::parallel_for(nz * ny, workers, [&](std::size_t row_begin, std::size_t row_end) {
        for (std::size_t row = row_begin; row < row_end; ++row) {
            const std::size_t k = row / ny;
            const std::size_t j = row % ny;
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t index = i + nx * (j + ny * k);
                const Scalar dhx_dz =
                    diff_backward(hx, index, sz, k, nz, axis_z.periodic, phase_z_inverse, inverse_dz);
                const Scalar dhz_dx =
                    diff_backward(hz, index, sx, i, nx, axis_x.periodic, phase_x_inverse, inverse_dx);
                Scalar stretched_z = dhx_dz;
                Scalar stretched_x = dhz_dx;
                if (!axis_z.periodic) {
                    const Scalar previous = psi_z[index];
                    stretched_z = dhx_dz * axis_z.inv_kappa_int[k] + previous;
                    psi_z[index] = axis_z.b_e_int[k] * previous + axis_z.c_e_int[k] * dhx_dz;
                }
                if (!axis_x.periodic) {
                    const Scalar previous = psi_x[index];
                    stretched_x = dhz_dx * axis_x.inv_kappa_int[i] + previous;
                    psi_x[index] = axis_x.b_e_int[i] * previous + axis_x.c_e_int[i] * dhz_dx;
                }
                const Scalar curl = stretched_z - stretched_x;
                if (dispersion) {
                    const Scalar previous = dispersion->jy[index];
                    const Scalar next = drude_a[index] * previous + drude_b[index] * ey[index];
                    ey[index] = ca[index] * ey[index] + cb[index] * (curl - 0.5 * (next + previous));
                    dispersion->jy[index] = next;
                } else {
                    ey[index] = ca[index] * ey[index] + cb[index] * curl;
                }
            }
        }
    });
}

template<class Scalar>
result<void> update_ez(const grid& domain, field_state<Scalar>& f, const detail::cpml_set& cpml,
               detail::cpml_memory<Scalar>& memory, const material_coefficients& coefficients,
               dispersion_state<Scalar>* dispersion, complex phase_x_inverse, complex phase_y_inverse,
               std::size_t workers) {
    const std::size_t nx = domain.cells[0], ny = domain.cells[1], nz = domain.cells[2];
    const std::size_t sx = 1, sy = nx;
    const double inverse_dx = 1.0 / domain.spacing[0], inverse_dy = 1.0 / domain.spacing[1];
    const auto& axis_x = cpml.axes[0];
    const auto& axis_y = cpml.axes[1];
    Scalar* ez = f.ez.data();
    const Scalar* hx = f.hx.data();
    const Scalar* hy = f.hy.data();
    Scalar* psi_x = memory.ez_x.empty() ? nullptr : memory.ez_x.data();
    Scalar* psi_y = memory.ez_y.empty() ? nullptr : memory.ez_y.data();
    const auto* ca = coefficients.ca[2].data();
    const auto* cb = coefficients.cb[2].data();
    const auto* drude_a = coefficients.has_drude ? coefficients.drude_a[2].data() : nullptr;
    const auto* drude_b = coefficients.has_drude ? coefficients.drude_b[2].data() : nullptr;

    return detail::parallel_for(nz * ny, workers, [&](std::size_t row_begin, std::size_t row_end) {
        for (std::size_t row = row_begin; row < row_end; ++row) {
            const std::size_t k = row / ny;
            const std::size_t j = row % ny;
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t index = i + nx * (j + ny * k);
                const Scalar dhy_dx =
                    diff_backward(hy, index, sx, i, nx, axis_x.periodic, phase_x_inverse, inverse_dx);
                const Scalar dhx_dy =
                    diff_backward(hx, index, sy, j, ny, axis_y.periodic, phase_y_inverse, inverse_dy);
                Scalar stretched_x = dhy_dx;
                Scalar stretched_y = dhx_dy;
                if (!axis_x.periodic) {
                    const Scalar previous = psi_x[index];
                    stretched_x = dhy_dx * axis_x.inv_kappa_int[i] + previous;
                    psi_x[index] = axis_x.b_e_int[i] * previous + axis_x.c_e_int[i] * dhy_dx;
                }
                if (!axis_y.periodic) {
                    const Scalar previous = psi_y[index];
                    stretched_y = dhx_dy * axis_y.inv_kappa_int[j] + previous;
                    psi_y[index] = axis_y.b_e_int[j] * previous + axis_y.c_e_int[j] * dhx_dy;
                }
                const Scalar curl = stretched_x - stretched_y;
                if (dispersion) {
                    const Scalar previous = dispersion->jz[index];
                    const Scalar next = drude_a[index] * previous + drude_b[index] * ez[index];
                    ez[index] = ca[index] * ez[index] + cb[index] * (curl - 0.5 * (next + previous));
                    dispersion->jz[index] = next;
                } else {
                    ez[index] = ca[index] * ez[index] + cb[index] * curl;
                }
            }
        }
    });
}

template<class Scalar>
result<void> update_hx(const grid& domain, field_state<Scalar>& f, const detail::cpml_set& cpml,
               detail::cpml_memory<Scalar>& memory, const material_coefficients& coefficients,
               complex phase_y, complex phase_z, std::size_t workers) {
    const std::size_t nx = domain.cells[0], ny = domain.cells[1], nz = domain.cells[2];
    const std::size_t sy = nx, sz = nx * ny;
    const double inverse_dy = 1.0 / domain.spacing[1], inverse_dz = 1.0 / domain.spacing[2];
    const auto& axis_y = cpml.axes[1];
    const auto& axis_z = cpml.axes[2];
    Scalar* hx = f.hx.data();
    const Scalar* ey = f.ey.data();
    const Scalar* ez = f.ez.data();
    Scalar* psi_y = memory.hx_y.empty() ? nullptr : memory.hx_y.data();
    Scalar* psi_z = memory.hx_z.empty() ? nullptr : memory.hx_z.data();
    const auto* cb = coefficients.h_cb[0].data();

    return detail::parallel_for(nz * ny, workers, [&](std::size_t row_begin, std::size_t row_end) {
        for (std::size_t row = row_begin; row < row_end; ++row) {
            const std::size_t k = row / ny;
            const std::size_t j = row % ny;
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t index = i + nx * (j + ny * k);
                const Scalar dey_dz =
                    diff_forward(ey, index, sz, k, nz, axis_z.periodic, phase_z, inverse_dz);
                const Scalar dez_dy =
                    diff_forward(ez, index, sy, j, ny, axis_y.periodic, phase_y, inverse_dy);
                Scalar stretched_z = dey_dz;
                Scalar stretched_y = dez_dy;
                if (!axis_z.periodic) {
                    const Scalar previous = psi_z[index];
                    stretched_z = dey_dz * axis_z.inv_kappa_half[k] + previous;
                    psi_z[index] = axis_z.b_h_half[k] * previous + axis_z.c_h_half[k] * dey_dz;
                }
                if (!axis_y.periodic) {
                    const Scalar previous = psi_y[index];
                    stretched_y = dez_dy * axis_y.inv_kappa_half[j] + previous;
                    psi_y[index] = axis_y.b_h_half[j] * previous + axis_y.c_h_half[j] * dez_dy;
                }
                hx[index] += cb[index] * (stretched_z - stretched_y);
            }
        }
    });
}

template<class Scalar>
result<void> update_hy(const grid& domain, field_state<Scalar>& f, const detail::cpml_set& cpml,
               detail::cpml_memory<Scalar>& memory, const material_coefficients& coefficients,
               complex phase_x, complex phase_z, std::size_t workers) {
    const std::size_t nx = domain.cells[0], ny = domain.cells[1], nz = domain.cells[2];
    const std::size_t sx = 1, sz = nx * ny;
    const double inverse_dx = 1.0 / domain.spacing[0], inverse_dz = 1.0 / domain.spacing[2];
    const auto& axis_x = cpml.axes[0];
    const auto& axis_z = cpml.axes[2];
    Scalar* hy = f.hy.data();
    const Scalar* ez = f.ez.data();
    const Scalar* ex = f.ex.data();
    Scalar* psi_x = memory.hy_x.empty() ? nullptr : memory.hy_x.data();
    Scalar* psi_z = memory.hy_z.empty() ? nullptr : memory.hy_z.data();
    const auto* cb = coefficients.h_cb[1].data();

    return detail::parallel_for(nz * ny, workers, [&](std::size_t row_begin, std::size_t row_end) {
        for (std::size_t row = row_begin; row < row_end; ++row) {
            const std::size_t k = row / ny;
            const std::size_t j = row % ny;
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t index = i + nx * (j + ny * k);
                const Scalar dez_dx =
                    diff_forward(ez, index, sx, i, nx, axis_x.periodic, phase_x, inverse_dx);
                const Scalar dex_dz =
                    diff_forward(ex, index, sz, k, nz, axis_z.periodic, phase_z, inverse_dz);
                Scalar stretched_x = dez_dx;
                Scalar stretched_z = dex_dz;
                if (!axis_x.periodic) {
                    const Scalar previous = psi_x[index];
                    stretched_x = dez_dx * axis_x.inv_kappa_half[i] + previous;
                    psi_x[index] = axis_x.b_h_half[i] * previous + axis_x.c_h_half[i] * dez_dx;
                }
                if (!axis_z.periodic) {
                    const Scalar previous = psi_z[index];
                    stretched_z = dex_dz * axis_z.inv_kappa_half[k] + previous;
                    psi_z[index] = axis_z.b_h_half[k] * previous + axis_z.c_h_half[k] * dex_dz;
                }
                hy[index] += cb[index] * (stretched_x - stretched_z);
            }
        }
    });
}

template<class Scalar>
result<void> update_hz(const grid& domain, field_state<Scalar>& f, const detail::cpml_set& cpml,
               detail::cpml_memory<Scalar>& memory, const material_coefficients& coefficients,
               complex phase_x, complex phase_y, std::size_t workers) {
    const std::size_t nx = domain.cells[0], ny = domain.cells[1], nz = domain.cells[2];
    const std::size_t sx = 1, sy = nx;
    const double inverse_dx = 1.0 / domain.spacing[0], inverse_dy = 1.0 / domain.spacing[1];
    const auto& axis_x = cpml.axes[0];
    const auto& axis_y = cpml.axes[1];
    Scalar* hz = f.hz.data();
    const Scalar* ex = f.ex.data();
    const Scalar* ey = f.ey.data();
    Scalar* psi_x = memory.hz_x.empty() ? nullptr : memory.hz_x.data();
    Scalar* psi_y = memory.hz_y.empty() ? nullptr : memory.hz_y.data();
    const auto* cb = coefficients.h_cb[2].data();

    return detail::parallel_for(nz * ny, workers, [&](std::size_t row_begin, std::size_t row_end) {
        for (std::size_t row = row_begin; row < row_end; ++row) {
            const std::size_t k = row / ny;
            const std::size_t j = row % ny;
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t index = i + nx * (j + ny * k);
                const Scalar dex_dy =
                    diff_forward(ex, index, sy, j, ny, axis_y.periodic, phase_y, inverse_dy);
                const Scalar dey_dx =
                    diff_forward(ey, index, sx, i, nx, axis_x.periodic, phase_x, inverse_dx);
                Scalar stretched_y = dex_dy;
                Scalar stretched_x = dey_dx;
                if (!axis_y.periodic) {
                    const Scalar previous = psi_y[index];
                    stretched_y = dex_dy * axis_y.inv_kappa_half[j] + previous;
                    psi_y[index] = axis_y.b_h_half[j] * previous + axis_y.c_h_half[j] * dex_dy;
                }
                if (!axis_x.periodic) {
                    const Scalar previous = psi_x[index];
                    stretched_x = dey_dx * axis_x.inv_kappa_half[i] + previous;
                    psi_x[index] = axis_x.b_h_half[i] * previous + axis_x.c_h_half[i] * dey_dx;
                }
                hz[index] += cb[index] * (stretched_y - stretched_x);
            }
        }
    });
}

template<class Scalar>
Scalar carrier_value(double phase) noexcept {
    if constexpr (std::is_same_v<Scalar, double>)
        return std::cos(phase);
    else
        return Scalar(std::cos(phase), std::sin(phase));
}

// Initialize a pulsed plane wave in the whole domain. E is placed at t = 0 and
// H at t = -dt/2, matching the leapfrog start of the time loop. The pulse is a
// one-way wave: no backward wave and no source-plane current discontinuity.
template<class Scalar>
result<void> initialize_fields(field_state<Scalar>& fields, const source_info& source, const grid& domain,
                               double dt, std::size_t workers) {
    const std::size_t nx = domain.cells[0], ny = domain.cells[1], nz = domain.cells[2];
    const std::array<Scalar*, 6> data{fields.ex.data(), fields.ey.data(), fields.ez.data(),
                                      fields.hx.data(), fields.hy.data(), fields.hz.data()};
    const double dx = domain.spacing[0], dy = domain.spacing[1], dz = domain.spacing[2];
    const double x0 = domain.origin[0], y0 = domain.origin[1], z0 = domain.origin[2];
    const double half_time_phase = 0.5 * source.omega * dt;
    return detail::parallel_for(nz * ny, workers, [&](std::size_t row_begin, std::size_t row_end) {
        for (std::size_t row = row_begin; row < row_end; ++row) {
            const std::size_t k = row / ny;
            const std::size_t j = row % ny;
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t index = i + nx * (j + ny * k);
                for (std::size_t c = 0; c < 6; ++c) {
                    const auto component_value = static_cast<component>(c);
                    const auto offset = component_offset(component_value);
                    const double x = x0 + (static_cast<double>(i) + offset[0]) * dx;
                    const double y = y0 + (static_cast<double>(j) + offset[1]) * dy;
                    const double z = z0 + (static_cast<double>(k) + offset[2]) * dz;
                    const double normalized = (z - source.z0) / source.sigma_z;
                    if (std::abs(normalized) > 7.0) {
                        data[c][index] = Scalar{};
                        continue;
                    }
                    const double envelope = std::exp(-0.5 * normalized * normalized);
                    double phase = source.kx * x + source.ky * y + source.kz * z + source.phase_offset;
                    if (c >= 3) phase += half_time_phase;
                    if (c < 3) {
                        data[c][index] = source.e0[c] * envelope * carrier_value<Scalar>(phase);
                    } else {
                        // Include the envelope-derivative term so the initial
                        // condition satisfies Maxwell's equations for a pulsed
                        // plane wave, not just for a monochromatic one.
                        const double slope = -(z - source.z0) / (source.sigma_z * source.sigma_z);
                        if constexpr (std::is_same_v<Scalar, double>) {
                            data[c][index] = envelope * (source.h0[c - 3] * std::cos(phase) +
                                                         slope * source.h1[c - 3] * std::sin(phase));
                        } else {
                            const complex correction{0.0, -slope * source.h1[c - 3]};
                            data[c][index] = envelope * (source.h0[c - 3] + correction) *
                                             carrier_value<Scalar>(phase);
                        }
                    }
                }
            }
        }
    });
}

template<class Scalar>
void accumulate_monitors(const grid& domain, const std::vector<monitor>& monitors,
                         std::vector<monitor_result>& accumulators,
                         const field_state<Scalar>& fields, double time_electric, double time_magnetic,
                         bool complex_mode) {
    const std::size_t nx = domain.cells[0], ny = domain.cells[1], nz = domain.cells[2];
    const std::array<const Scalar*, 6> data{fields.ex.data(), fields.ey.data(), fields.ez.data(),
                                            fields.hx.data(), fields.hy.data(), fields.hz.data()};
    for (std::size_t monitor_index = 0; monitor_index < monitors.size(); ++monitor_index) {
        auto& target = accumulators[monitor_index];
        const auto& definition = monitors[monitor_index];
        bool electric_inside = true;
        bool magnetic_inside = true;
        if (definition.time_window) {
            const auto& window = *definition.time_window;
            electric_inside = time_electric >= window[0] && time_electric <= window[1];
            magnetic_inside = time_magnetic >= window[0] && time_magnetic <= window[1];
        }
        if (!electric_inside && !magnetic_inside) continue;
        const std::size_t cells = target.cell_count();
        const std::size_t frequency_count = target.frequencies.size();
        std::vector<complex> phase_electric(frequency_count);
        std::vector<complex> phase_magnetic(frequency_count);
        for (std::size_t f = 0; f < frequency_count; ++f) {
            const double omega = 2.0 * std::numbers::pi * target.frequencies[f];
            if (complex_mode) {
                phase_electric[f] = std::exp(complex{0.0, omega * time_electric});
                phase_magnetic[f] = std::exp(complex{0.0, omega * time_magnetic});
            } else {
                phase_electric[f] = 2.0 * std::exp(complex{0.0, -omega * time_electric});
                phase_magnetic[f] = 2.0 * std::exp(complex{0.0, -omega * time_magnetic});
            }
        }

        auto accumulate_cell = [&](std::size_t index, std::size_t plane_index) {
            for (std::size_t f = 0; f < frequency_count; ++f) {
                for (std::size_t c = 0; c < 6; ++c) {
                    if (c < 3 ? !electric_inside : !magnetic_inside) continue;
                    const complex phase = c < 3 ? phase_electric[f] : phase_magnetic[f];
                    auto& slot = target.values[(f * 6 + c) * cells + plane_index];
                    slot += complex(data[c][index]) * phase;
                }
            }
        };

        switch (definition.normal) {
        case axis::z:
            for (std::size_t j = 0; j < ny; ++j)
                for (std::size_t i = 0; i < nx; ++i)
                    accumulate_cell(i + nx * (j + ny * definition.index), i + nx * j);
            break;
        case axis::x:
            for (std::size_t k = 0; k < nz; ++k)
                for (std::size_t j = 0; j < ny; ++j)
                    accumulate_cell(definition.index + nx * (j + ny * k), j + ny * k);
            break;
        case axis::y:
            for (std::size_t k = 0; k < nz; ++k)
                for (std::size_t i = 0; i < nx; ++i)
                    accumulate_cell(i + nx * (definition.index + ny * k), i + nx * k);
            break;
        }
    }
}

result<material_coefficients> build_coefficients(const problem& value, const material_map& map, double dt) noexcept {
    return detail::guard([&]() -> result<material_coefficients> {
        material_coefficients out;
        bool has_drude = false;
        for (const auto& entry : value.materials.all())
            if (entry.dispersive()) has_drude = true;
        out.has_drude = has_drude;

        const std::size_t count = value.grid.size();
        for (std::size_t c = 0; c < 3; ++c) {
            out.ca[c].resize(count);
            out.cb[c].resize(count);
            if (has_drude) {
                out.drude_a[c].assign(count, 0.0);
                out.drude_b[c].assign(count, 0.0);
            }
            const auto ids = map[static_cast<component>(c)];
            for (std::size_t index = 0; index < count; ++index) {
                const material* entry = value.materials.at(ids[index]);
                if (!entry) return fail(errc::missing_material);
                const double epsilon = entry->epsilon_inf;
                const double loss = entry->conductivity * dt / (2.0 * vacuum_permittivity * epsilon);
                out.ca[c][index] = (1.0 - loss) / (1.0 + loss);
                out.cb[c][index] = (dt / (vacuum_permittivity * epsilon)) / (1.0 + loss);
                if (entry->drude) {
                    const double gamma = entry->drude->collision_rate;
                    const double denominator = 1.0 + 0.5 * gamma * dt;
                    out.drude_a[c][index] = (1.0 - 0.5 * gamma * dt) / denominator;
                    out.drude_b[c][index] = dt * vacuum_permittivity * entry->drude->plasma_frequency *
                                            entry->drude->plasma_frequency / denominator;
                }
            }
        }
        for (std::size_t c = 0; c < 3; ++c) {
            out.h_cb[c].resize(count);
            const auto ids = map[static_cast<component>(3 + c)];
            for (std::size_t index = 0; index < count; ++index) {
                const material* entry = value.materials.at(ids[index]);
                if (!entry) return fail(errc::missing_material);
                out.h_cb[c][index] = dt / (vacuum_permeability * entry->mu_r);
            }
        }
        return out;
    });
}

result<source_info> make_source_info(const problem& value, const material_map& map) noexcept {
    return detail::guard([&]() -> result<source_info> {
        source_info out;
        const auto& domain = value.grid;
        const auto& source = value.source;
        const std::size_t nx = domain.cells[0], ny = domain.cells[1], nz = domain.cells[2];
        const double dz = domain.spacing[2];
        out.z0 = source.z;
        out.pulse_sigma = pulse_sigma(source);

        const auto center_k = static_cast<std::size_t>(
            std::clamp<long long>(static_cast<long long>((source.z - domain.origin[2]) / dz), 0,
                                  static_cast<long long>(nz) - 1));
        const std::uint32_t id = map.at(component::ex, nx * ny * center_k);
        const material* medium = value.materials.at(id);
        if (!medium) return fail(errc::missing_material);
        if (medium->dispersive() || medium->conductivity != 0) return fail(errc::unsupported);

        const double k0 = 2.0 * std::numbers::pi / source.wavelength;
        const double index = std::sqrt(medium->epsilon_inf * medium->mu_r);
        out.sigma_z = speed_of_light / index * std::cos(source.theta) * out.pulse_sigma;
        if (!(out.sigma_z > 0)) return fail(errc::invalid_argument);

        const double z_low = source.z - 4.0 * out.sigma_z;
        const double z_high = source.z + 4.0 * out.sigma_z;
        const double domain_low = domain.origin[2];
        const double domain_high = domain.origin[2] + static_cast<double>(nz) * dz;
        const double direction_sign = source.from_below ? 1.0 : -1.0;
        const double forward_edge = source.z + direction_sign * 4.0 * out.sigma_z;
        if (source.z < domain_low || source.z > domain_high) return fail(errc::out_of_range);
        if (forward_edge < domain_low || forward_edge > domain_high) return fail(errc::out_of_range);
        if (!is_periodic(value.boundaries, 2)) {
            const double pml_low = domain_low + static_cast<double>(value.boundaries.cpml.layers) * dz;
            const double pml_high = domain_high - static_cast<double>(value.boundaries.cpml.layers) * dz;
            if (source.z < pml_low || source.z > pml_high) return fail(errc::source_inside_pml);
            if (direction_sign > 0 ? forward_edge > pml_high : forward_edge < pml_low)
                return fail(errc::source_inside_pml);
        }

        // The pulse must live in one uniform, lossless material.
        const auto first_k = static_cast<std::size_t>(
            std::clamp<long long>(static_cast<long long>(std::floor((z_low - domain.origin[2]) / dz)), 0,
                                  static_cast<long long>(nz) - 1));
        const auto last_k = static_cast<std::size_t>(
            std::clamp<long long>(static_cast<long long>(std::ceil((z_high - domain.origin[2]) / dz)), 0,
                                  static_cast<long long>(nz) - 1));
        for (std::size_t k = first_k; k <= last_k; ++k) {
            for (std::size_t j = 0; j < ny; ++j) {
                for (std::size_t i = 0; i < nx; ++i) {
                    const std::size_t index = i + nx * (j + ny * k);
                    for (std::size_t c = 0; c < 3; ++c)
                        if (map.at(static_cast<component>(c), index) != id)
                            return fail(errc::source_plane_not_uniform);
                }
            }
        }

        const double direction = source.from_below ? 1.0 : -1.0;
        const double transverse = index * k0 * std::sin(source.theta);
        out.kx = transverse * std::cos(source.phi);
        out.ky = transverse * std::sin(source.phi);
        out.kz = direction * index * k0 * std::cos(source.theta);
        out.k_magnitude = index * k0;
        out.omega = 2.0 * std::numbers::pi * frequency(source);
        out.frequency = frequency(source);
        out.phase_offset = source.phase;

        const common::vector3d s{-std::sin(source.phi), std::cos(source.phi), 0.0};
        const common::vector3d k_hat{out.kx / out.k_magnitude, out.ky / out.k_magnitude,
                                     out.kz / out.k_magnitude};
        const common::vector3d p = common::cross(s, k_hat);
        const common::vector3d e0 = source.amplitude * (std::cos(source.polarization) * s +
                                                        std::sin(source.polarization) * p);
        const common::vector3d k_vector{out.kx, out.ky, out.kz};
        const common::vector3d h0 = common::cross(k_vector, e0) /
                                    (out.omega * vacuum_permeability * medium->mu_r);
        const common::vector3d h1 = common::cross(common::vector3d{0.0, 0.0, 1.0}, e0) /
                                    (out.omega * vacuum_permeability * medium->mu_r);
        out.e0 = {e0[0], e0[1], e0[2]};
        out.h0 = {h0[0], h0[1], h0[2]};
        out.h1 = {h1[0], h1[1], h1[2]};

        out.complex_mode = false;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (!is_periodic(value.boundaries, axis)) continue;
            const double k = axis == 0 ? out.kx : axis == 1 ? out.ky : out.kz;
            if (std::abs(k) > 1e-12 * out.k_magnitude) out.complex_mode = true;
        }
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (out.complex_mode && is_periodic(value.boundaries, axis)) {
                const double k = axis == 0 ? out.kx : axis == 1 ? out.ky : out.kz;
                const double length = static_cast<double>(domain.cells[axis]) * domain.spacing[axis];
                out.phase[axis] = std::exp(complex{0.0, k * length});
                out.phase_inv[axis] = complex{1.0, 0.0} / out.phase[axis];
            }
        }
        return out;
    });
}

template<class Scalar>
result<run_result> run_impl(const problem& value, const source_info& source,
                            const detail::cpml_set& cpml, const material_coefficients& coefficients,
                            double dt, std::size_t steps) {
    return detail::guard([&]() -> result<run_result> {
        const auto& domain = value.grid;
        const std::size_t count = domain.size();
        field_state<Scalar> fields;
        allocate(fields, count);
        detail::cpml_memory<Scalar> memory;
        detail::allocate(memory, cpml, count);

        dispersion_state<Scalar> dispersion;
        if (coefficients.has_drude) {
            dispersion.jx.assign(count, Scalar{});
            dispersion.jy.assign(count, Scalar{});
            dispersion.jz.assign(count, Scalar{});
        }
        auto initialized = initialize_fields(fields, source, domain, dt, value.run.workers);
        if (!initialized) return fail(initialized.error());

        std::vector<monitor_result> accumulators;
        accumulators.reserve(value.monitors.size());
        for (const auto& definition : value.monitors) {
            monitor_result result;
            result.name = definition.name;
            result.normal = definition.normal;
            result.index = definition.index;
            result.shape = plane_shape(domain, definition.normal);
            result.spacing = plane_spacing(domain, definition.normal);
            result.origin = plane_origin(domain, definition.normal);
            result.normal_spacing = domain.spacing[static_cast<std::size_t>(definition.normal)];
            result.frequencies = definition.frequencies.empty()
                                     ? std::vector<double>{source.frequency}
                                     : definition.frequencies;
            const std::size_t cells = result.cell_count();
            result.values.assign(result.frequencies.size() * 6 * cells, complex{});
            accumulators.push_back(std::move(result));
        }

        const std::size_t workers = value.run.workers;
        const auto start = std::chrono::steady_clock::now();
        for (std::size_t step = 0; step < steps; ++step) {
            const double time_magnetic = (static_cast<double>(step) + 0.5) * dt;
            const double time_electric = (static_cast<double>(step) + 1.0) * dt;
            if (auto status = update_hx(domain, fields, cpml, memory, coefficients, source.phase[1],
                                        source.phase[2], workers);
                !status)
                return fail(status.error());
            if (auto status = update_hy(domain, fields, cpml, memory, coefficients, source.phase[0],
                                        source.phase[2], workers);
                !status)
                return fail(status.error());
            if (auto status = update_hz(domain, fields, cpml, memory, coefficients, source.phase[0],
                                        source.phase[1], workers);
                !status)
                return fail(status.error());
            if (auto status = update_ex(domain, fields, cpml, memory, coefficients,
                                        coefficients.has_drude ? &dispersion : nullptr, source.phase_inv[1],
                                        source.phase_inv[2], workers);
                !status)
                return fail(status.error());
            if (auto status = update_ey(domain, fields, cpml, memory, coefficients,
                                        coefficients.has_drude ? &dispersion : nullptr, source.phase_inv[0],
                                        source.phase_inv[2], workers);
                !status)
                return fail(status.error());
            if (auto status = update_ez(domain, fields, cpml, memory, coefficients,
                                        coefficients.has_drude ? &dispersion : nullptr, source.phase_inv[0],
                                        source.phase_inv[1], workers);
                !status)
                return fail(status.error());
            accumulate_monitors(domain, value.monitors, accumulators, fields, time_electric,
                                time_magnetic, source.complex_mode);
        }
        const auto finish = std::chrono::steady_clock::now();

        // The DFT sum of a Gaussian pulse is its envelope integral. Dividing
        // by the analytic integral gives unit spectrum coefficients for a unit
        // amplitude plane wave, independent of the monitor position.
        const double envelope_integral = std::sqrt(2.0 * std::numbers::pi) * source.pulse_sigma / dt;
        if (envelope_integral > 0) {
            for (auto& target : accumulators) {
                const std::size_t cells = target.cell_count();
                const std::size_t frequency_count = target.frequencies.size();
                for (std::size_t f = 0; f < frequency_count; ++f)
                    for (std::size_t c = 0; c < 6; ++c) {
                        auto* data = target.values.data() + (f * 6 + c) * cells;
                        for (std::size_t cell = 0; cell < cells; ++cell) data[cell] /= envelope_integral;
                    }
            }
        }

        run_result result;
        result.dt = dt;
        result.steps = steps;
        result.frequency = source.frequency;
        result.source_wavevector = {source.kx, source.ky, source.kz};
        result.workers = detail::worker_count(workers);
        result.monitors = std::move(accumulators);
        result.elapsed_seconds = std::chrono::duration<double>(finish - start).count();

        for (std::size_t index = 0; index < result.monitors.size(); ++index) {
            if (!value.monitors[index].spectrum) continue;
            spectrum_options options = *value.monitors[index].spectrum;
            if (!options.incident_wavevector) options.incident_wavevector = std::array<double, 2>{source.kx, source.ky};
            auto spectrum = compute_spectrum(result.monitors[index], options);
            if (!spectrum) return fail(spectrum.error());
            result.monitors[index].spectrum = std::move(*spectrum);
        }
        return result;
    });
}

} // namespace

result<run_result> solve(const problem& value) noexcept {
    return detail::guard([&]() -> result<run_result> {
        if (auto valid = validate(value.grid); !valid) return fail(valid.error());
        if (value.materials.size() == 0) return fail(errc::missing_material);
        if (auto valid = validate(value.geometry, value.materials); !valid) return fail(valid.error());
        if (auto valid = validate(value.boundaries, value.grid); !valid) return fail(valid.error());
        if (auto valid = validate(value.source); !valid) return fail(valid.error());
        if (auto valid = validate(value.run); !valid) return fail(valid.error());

        for (const auto& monitor : value.monitors) {
            if (auto valid = validate(monitor, value.grid); !valid) return fail(valid.error());
            const auto normal = static_cast<std::size_t>(monitor.normal);
            if (in_cpml(value.boundaries, normal, monitor.index, value.grid.cells[normal]))
                return fail(errc::monitor_inside_pml);
        }

        auto map = voxelize(value.geometry, value.materials, value.grid);
        if (!map) return fail(map.error());
        auto source = make_source_info(value, *map);
        if (!source) return fail(source.error());

        double minimum_epsilon = std::numeric_limits<double>::infinity();
        double minimum_mu = std::numeric_limits<double>::infinity();
        for (const auto& entry : value.materials.all()) {
            minimum_epsilon = std::min(minimum_epsilon, entry.epsilon_inf);
            minimum_mu = std::min(minimum_mu, entry.mu_r);
        }
        if (!(minimum_epsilon > 0) || !(minimum_mu > 0)) return fail(errc::invalid_configuration);
        const double maximum_speed = speed_of_light / std::sqrt(minimum_epsilon * minimum_mu);
        const auto& spacing = value.grid.spacing;
        const double inverse_length =
            std::sqrt(1.0 / (spacing[0] * spacing[0]) + 1.0 / (spacing[1] * spacing[1]) +
                      1.0 / (spacing[2] * spacing[2]));
        const double dt = value.run.courant / (maximum_speed * inverse_length);
        if (!std::isfinite(dt) || !(dt > 0)) return fail(errc::unstable);
        for (const auto& entry : value.materials.all())
            if (entry.drude && entry.drude->plasma_frequency * dt > 1.0) return fail(errc::unstable);

        std::size_t steps = 0;
        if (value.run.steps) {
            steps = *value.run.steps;
        } else {
            const double duration = value.run.duration * period(value.source);
            const double count = std::ceil(duration / dt);
            if (!std::isfinite(count) || count > static_cast<double>(std::numeric_limits<std::size_t>::max()))
                return fail(errc::out_of_range);
            steps = std::max<std::size_t>(1, static_cast<std::size_t>(count));
        }

        auto cpml = detail::make_cpml_set(value.grid, value.boundaries, dt);
        if (!cpml) return fail(cpml.error());
        auto coefficients = build_coefficients(value, *map, dt);
        if (!coefficients) return fail(coefficients.error());

        if (source->complex_mode)
            return run_impl<complex>(value, *source, *cpml, *coefficients, dt, steps);
        return run_impl<double>(value, *source, *cpml, *coefficients, dt, steps);
    });
}

} // namespace litho::fdtd
