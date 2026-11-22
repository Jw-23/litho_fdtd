#pragma once

#include <litho_fdtd/boundary.hpp>
#include <litho_fdtd/constants.hpp>
#include <litho_fdtd/error.hpp>
#include <litho_fdtd/grid.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace litho::fdtd::detail {

// Per-axis CPML coefficients at integer and half-integer node positions. The
// two lattices cover the staggered field locations. Periodic axes carry no
// coefficients and are handled by phase-shifted wrapping in the update loops.
struct cpml_axis {
    bool periodic = false;
    std::size_t count = 0;
    double spacing = 0;
    std::vector<double> inv_kappa_int, b_e_int, c_e_int, b_h_int, c_h_int;
    std::vector<double> inv_kappa_half, b_e_half, c_e_half, b_h_half, c_h_half;
    std::vector<std::uint8_t> pml_int, pml_half;
};

struct cpml_set {
    std::array<cpml_axis, 3> axes;
};

namespace cpml_detail {

struct profile {
    double sigma = 0;
    double kappa = 1;
    double alpha = 0;
    bool active = false;
};

inline profile sample(double coordinate, double length, double thickness, const cpml_options& options,
                      double sigma_max) noexcept {
    double depth = 0;
    if (coordinate < thickness)
        depth = thickness - coordinate;
    else if (coordinate > length - thickness)
        depth = coordinate - (length - thickness);
    const double u = std::clamp(depth / thickness, 0.0, 1.0);
    const double graded = std::pow(u, options.order);
    profile out;
    out.sigma = sigma_max * graded;
    out.kappa = 1.0 + (options.kappa_max - 1.0) * graded;
    out.alpha = options.alpha_max * sigma_max * (1.0 - u);
    out.active = depth > 0;
    return out;
}

struct coefficients {
    double inv_kappa;
    double b_e, c_e;
    double b_h, c_h;
};

inline coefficients make_coefficients(const profile& p, double dt) noexcept {
    const double sigma_over_kappa = p.sigma / p.kappa;
    const double b_e = std::exp(-(sigma_over_kappa + p.alpha) * dt / vacuum_permittivity);
    const double b_h = std::exp(-(sigma_over_kappa + p.alpha) * dt / vacuum_permeability);
    const double denominator = p.sigma * p.kappa + p.kappa * p.kappa * p.alpha;
    const double c_e = denominator == 0 ? 0.0 : p.sigma / denominator * (b_e - 1.0);
    const double c_h = denominator == 0 ? 0.0 : p.sigma / denominator * (b_h - 1.0);
    return {1.0 / p.kappa, b_e, c_e, b_h, c_h};
}

} // namespace cpml_detail

[[nodiscard]] inline result<cpml_axis> make_cpml_axis(const grid& domain, const boundaries& bc,
                                                      std::size_t axis, double dt) noexcept {
    cpml_axis out;
    out.periodic = is_periodic(bc, axis);
    out.count = domain.cells[axis];
    out.spacing = domain.spacing[axis];
    if (out.periodic) return out;

    const double length = static_cast<double>(out.count) * out.spacing;
    const double thickness = static_cast<double>(bc.cpml.layers) * out.spacing;
    if (!(thickness > 0) || 2.0 * thickness >= length) return fail(errc::invalid_configuration);
    const double sigma_max = -(bc.cpml.order + 1.0) * std::log(bc.cpml.reflection) /
                             (2.0 * vacuum_impedance * thickness);
    if (!std::isfinite(sigma_max) || !(sigma_max > 0)) return fail(errc::invalid_configuration);

    const std::size_t n = out.count;
    out.inv_kappa_int.resize(n);
    out.b_e_int.resize(n);
    out.c_e_int.resize(n);
    out.b_h_int.resize(n);
    out.c_h_int.resize(n);
    out.inv_kappa_half.resize(n);
    out.b_e_half.resize(n);
    out.c_e_half.resize(n);
    out.b_h_half.resize(n);
    out.c_h_half.resize(n);
    out.pml_int.resize(n);
    out.pml_half.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        const auto node_profile = cpml_detail::sample(static_cast<double>(i) * out.spacing, length, thickness,
                                                      bc.cpml, sigma_max);
        const auto half_profile = cpml_detail::sample((static_cast<double>(i) + 0.5) * out.spacing, length,
                                                      thickness, bc.cpml, sigma_max);
        const auto node = cpml_detail::make_coefficients(node_profile, dt);
        const auto half = cpml_detail::make_coefficients(half_profile, dt);
        out.pml_int[i] = node_profile.active ? 1 : 0;
        out.pml_half[i] = half_profile.active ? 1 : 0;
        out.inv_kappa_int[i] = node.inv_kappa;
        out.b_e_int[i] = node.b_e;
        out.c_e_int[i] = node.c_e;
        out.b_h_int[i] = node.b_h;
        out.c_h_int[i] = node.c_h;
        out.inv_kappa_half[i] = half.inv_kappa;
        out.b_e_half[i] = half.b_e;
        out.c_e_half[i] = half.c_e;
        out.b_h_half[i] = half.b_h;
        out.c_h_half[i] = half.c_h;
    }
    return out;
}

[[nodiscard]] inline result<cpml_set> make_cpml_set(const grid& domain, const boundaries& bc,
                                                    double dt) noexcept {
    cpml_set out;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        auto value = make_cpml_axis(domain, bc, axis, dt);
        if (!value) return fail(value.error());
        out.axes[axis] = std::move(*value);
    }
    return out;
}

// CPML auxiliary variables. A vector is left empty when its axis is periodic,
// which makes the corresponding branch a no-op and saves memory.
template<class Scalar>
struct cpml_memory {
    std::vector<Scalar> ex_y, ex_z;
    std::vector<Scalar> ey_x, ey_z;
    std::vector<Scalar> ez_x, ez_y;
    std::vector<Scalar> hx_y, hx_z;
    std::vector<Scalar> hy_x, hy_z;
    std::vector<Scalar> hz_x, hz_y;
};

template<class Scalar>
void allocate(cpml_memory<Scalar>& memory, const cpml_set& set, std::size_t count) {
    if (!set.axes[0].periodic) {
        memory.ey_x.resize(count);
        memory.ez_x.resize(count);
        memory.hy_x.resize(count);
        memory.hz_x.resize(count);
    }
    if (!set.axes[1].periodic) {
        memory.ex_y.resize(count);
        memory.ez_y.resize(count);
        memory.hx_y.resize(count);
        memory.hz_y.resize(count);
    }
    if (!set.axes[2].periodic) {
        memory.ex_z.resize(count);
        memory.ey_z.resize(count);
        memory.hx_z.resize(count);
        memory.hy_z.resize(count);
    }
}

} // namespace litho::fdtd::detail
