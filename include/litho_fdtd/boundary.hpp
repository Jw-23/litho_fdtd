#pragma once

#include <litho_fdtd/error.hpp>
#include <litho_fdtd/grid.hpp>

#include <array>
#include <cmath>
#include <cstddef>

namespace litho::fdtd {

enum class boundary_kind { cpml, periodic };

struct cpml_options {
    std::size_t layers = 10;
    double reflection = 1e-6; // target normal-incidence reflection
    double order = 3.0;       // polynomial grading exponent
    double kappa_max = 5.0;   // coordinate stretching at the outer edge
    double alpha_max = 0.05;  // CFS alpha as a fraction of sigma_max
};

struct boundaries {
    std::array<boundary_kind, 3> axes{boundary_kind::cpml, boundary_kind::cpml, boundary_kind::cpml};
    cpml_options cpml{};
};

[[nodiscard]] inline bool is_periodic(const boundaries& value, std::size_t axis) noexcept {
    return value.axes[axis] == boundary_kind::periodic;
}

[[nodiscard]] inline bool in_cpml(const boundaries& value, std::size_t axis, std::size_t index,
                                  std::size_t count) noexcept {
    if (is_periodic(value, axis) || value.cpml.layers == 0) return false;
    return index < value.cpml.layers || index + value.cpml.layers >= count;
}

[[nodiscard]] inline result<void> validate(const boundaries& value, const grid& domain) noexcept {
    bool any_cpml = false;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (value.axes[axis] == boundary_kind::periodic) continue;
        any_cpml = true;
        if (value.cpml.layers == 0) return fail(errc::invalid_configuration);
        if (2 * value.cpml.layers >= domain.cells[axis]) return fail(errc::invalid_configuration);
    }
    if (!any_cpml) return {};
    if (!std::isfinite(value.cpml.reflection) || !(value.cpml.reflection > 0) || !(value.cpml.reflection < 1))
        return fail(errc::invalid_configuration);
    if (!std::isfinite(value.cpml.order) || !(value.cpml.order > 0))
        return fail(errc::invalid_configuration);
    if (!std::isfinite(value.cpml.kappa_max) || !(value.cpml.kappa_max >= 1))
        return fail(errc::invalid_configuration);
    if (!std::isfinite(value.cpml.alpha_max) || value.cpml.alpha_max < 0)
        return fail(errc::invalid_configuration);
    return {};
}

} // namespace litho::fdtd
