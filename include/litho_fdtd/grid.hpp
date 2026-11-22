#pragma once

#include <litho_fdtd/error.hpp>
#include <litho_geometry/common.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>

namespace litho::fdtd {

// Geometry types live in litho::geometry::common; this alias keeps the solver
// API readable without a second set of point/vector types.
namespace common = litho::geometry::common;

// Uniform Cartesian grid. Fields are stored x-fastest over `cells` Yee cells.
// The physical domain is [origin, origin + cells * spacing].
struct grid {
    common::point3d origin{};
    std::array<std::size_t, 3> cells{};
    std::array<double, 3> spacing{};

    [[nodiscard]] constexpr std::size_t size() const noexcept {
        return cells[0] * cells[1] * cells[2];
    }
    [[nodiscard]] constexpr std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
        return i + cells[0] * (j + cells[1] * k);
    }
    [[nodiscard]] constexpr std::array<double, 3> extent() const noexcept {
        return {cells[0] * spacing[0], cells[1] * spacing[1], cells[2] * spacing[2]};
    }
    [[nodiscard]] constexpr common::point3d cell_center(std::size_t i, std::size_t j, std::size_t k) const noexcept {
        return {origin[0] + (static_cast<double>(i) + 0.5) * spacing[0],
                origin[1] + (static_cast<double>(j) + 0.5) * spacing[1],
                origin[2] + (static_cast<double>(k) + 0.5) * spacing[2]};
    }
    [[nodiscard]] constexpr double coordinate(std::size_t axis, double index) const noexcept {
        return origin[axis] + index * spacing[axis];
    }
    [[nodiscard]] constexpr common::box3d bounds() const noexcept {
        const auto e = extent();
        return {{origin[0], origin[1], origin[2]},
                {origin[0] + e[0], origin[1] + e[1], origin[2] + e[2]}};
    }
};

[[nodiscard]] inline result<void> validate(const grid& value) noexcept {
    if (!common::is_finite(value.origin)) return fail(errc::invalid_argument);
    for (std::size_t axis = 0; axis < 3; ++axis) {
        if (value.cells[axis] == 0) return fail(errc::invalid_argument);
        if (!std::isfinite(value.spacing[axis]) || !(value.spacing[axis] > 0))
            return fail(errc::invalid_argument);
    }
    const double max = static_cast<double>(std::numeric_limits<std::size_t>::max());
    if (static_cast<double>(value.cells[0]) * static_cast<double>(value.cells[1])
            * static_cast<double>(value.cells[2]) > max * 0.5)
        return fail(errc::out_of_range);
    if (!std::isfinite(value.extent()[0]) || !std::isfinite(value.extent()[1]) || !std::isfinite(value.extent()[2]))
        return fail(errc::invalid_argument);
    return {};
}

enum class axis : std::size_t { x = 0, y = 1, z = 2 };

// Field components are staggered on the lower edges and faces of each Yee cell.
enum class component : std::size_t { ex = 0, ey = 1, ez = 2, hx = 3, hy = 4, hz = 5 };

[[nodiscard]] constexpr std::array<double, 3> component_offset(component value) noexcept {
    switch (value) {
    case component::ex: return {0.5, 0.0, 0.0};
    case component::ey: return {0.0, 0.5, 0.0};
    case component::ez: return {0.0, 0.0, 0.5};
    case component::hx: return {0.0, 0.5, 0.5};
    case component::hy: return {0.5, 0.0, 0.5};
    case component::hz: return {0.5, 0.5, 0.0};
    }
    return {};
}

[[nodiscard]] constexpr bool is_electric(component value) noexcept {
    return value == component::ex || value == component::ey || value == component::ez;
}

} // namespace litho::fdtd
