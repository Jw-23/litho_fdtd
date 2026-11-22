#pragma once

#include <litho_fdtd/grid.hpp>
#include <litho_fdtd/material.hpp>
#include <litho_geometry/common.hpp>

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace litho::fdtd {

// A planar layer covers its z interval over the full xy extent. Layers are
// applied in order; a later layer overrides an earlier one.
struct planar_layer {
    double z_min = 0;
    double z_max = 0;
    std::size_t material = 0;
};

// A patterned solid overrides the planar stack. Solids are applied in order;
// a later solid overrides an earlier one.
struct patterned_solid {
    common::solid3d shape;
    std::size_t material = 0;
};

struct structure {
    std::size_t background = 0;
    std::vector<planar_layer> layers;
    std::vector<patterned_solid> solids;
};

// Material id at every staggered field location. Memory is intentionally kept
// as ids rather than coefficients so setup and diagnostics share one map.
struct material_map {
    std::array<std::vector<std::uint32_t>, 6> ids;

    [[nodiscard]] std::span<const std::uint32_t> operator[](component value) const noexcept {
        return ids[static_cast<std::size_t>(value)];
    }
    [[nodiscard]] std::span<std::uint32_t> operator[](component value) noexcept {
        return ids[static_cast<std::size_t>(value)];
    }
    [[nodiscard]] std::uint32_t at(component value, std::size_t index) const noexcept {
        return ids[static_cast<std::size_t>(value)][index];
    }
};

[[nodiscard]] result<void> validate(const structure& value, const material_library& materials) noexcept;

// Rasterize the structure onto the six staggered field lattices. This is a
// slice-oriented operation: for each z plane the solid is reduced to a 2D
// polygon and tested against the field locations on that plane.
[[nodiscard]] result<material_map> voxelize(const structure& value,
                                            const material_library& materials,
                                            const grid& domain) noexcept;

} // namespace litho::fdtd
