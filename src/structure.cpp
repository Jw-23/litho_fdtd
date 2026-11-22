#include <litho_fdtd/structure.hpp>

#include <boost/geometry.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace litho::fdtd {

namespace {

struct index_range {
    std::size_t begin = 0;
    std::size_t end = 0;
    bool empty() const noexcept { return begin >= end; }
};

index_range axis_range(double minimum,
                       double maximum,
                       std::size_t count,
                       double origin,
                       double spacing,
                       double offset) noexcept {
    const double low = (minimum - origin) / spacing - offset;
    const double high = (maximum - origin) / spacing - offset;
    if (!std::isfinite(low) || !std::isfinite(high)) return {};
    const auto first = static_cast<long long>(std::ceil(low - 1e-9));
    const auto last = static_cast<long long>(std::floor(high + 1e-9));
    const auto last_index = static_cast<long long>(count) - 1;
    const auto begin = std::clamp<long long>(first, 0, last_index);
    const auto end = std::clamp<long long>(last + 1, 0, last_index + 1);
    if (begin >= end) return {};
    return {static_cast<std::size_t>(begin), static_cast<std::size_t>(end)};
}

} // namespace

result<void> validate(const structure& value, const material_library& materials) noexcept {
    if (materials.size() == 0) return fail(errc::missing_material);
    if (value.background >= materials.size()) return fail(errc::missing_material);
    for (const auto& entry : materials.all())
        if (auto valid = validate(entry); !valid) return fail(valid.error());

    for (const auto& layer : value.layers) {
        if (!std::isfinite(layer.z_min) || !std::isfinite(layer.z_max) || !(layer.z_min < layer.z_max))
            return fail(errc::invalid_geometry);
        if (layer.material >= materials.size()) return fail(errc::missing_material);
    }
    for (const auto& solid : value.solids) {
        if (solid.material >= materials.size()) return fail(errc::missing_material);
        if (auto bounds = common::bounds(solid.shape); !bounds) return fail(errc::invalid_geometry);
    }
    return {};
}

result<material_map> voxelize(const structure& value,
                              const material_library& materials,
                              const grid& domain) noexcept {
    return detail::guard([&]() -> result<material_map> {
        if (auto valid = validate(domain); !valid) return fail(valid.error());
        if (auto valid = validate(value, materials); !valid) return fail(valid.error());

        material_map map;
        const std::size_t count = domain.size();
        for (auto& ids : map.ids) ids.assign(count, static_cast<std::uint32_t>(value.background));

        const std::size_t nx = domain.cells[0];
        const std::size_t ny = domain.cells[1];
        const std::size_t nz = domain.cells[2];
        const double dz = domain.spacing[2];

        // Planar layers: one material id per z slice, identical over xy.
        for (std::size_t c = 0; c < 6; ++c) {
            const auto part = static_cast<component>(c);
            const double offset_z = component_offset(part)[2];
            auto& ids = map.ids[c];
            for (std::size_t k = 0; k < nz; ++k) {
                const double z = domain.origin[2] + (static_cast<double>(k) + offset_z) * dz;
                std::uint32_t id = static_cast<std::uint32_t>(value.background);
                for (const auto& layer : value.layers)
                    if (z >= layer.z_min && z < layer.z_max) id = static_cast<std::uint32_t>(layer.material);
                std::fill_n(ids.begin() + static_cast<std::ptrdiff_t>(k * nx * ny),
                            static_cast<std::ptrdiff_t>(nx * ny), id);
            }
        }

        // Patterned solids: slice the solid at each z plane and rasterize the
        // resulting polygon on that plane's staggered xy lattice.
        for (const auto& solid : value.solids) {
            const auto bounds = common::bounds(solid.shape);
            if (!bounds) return fail(errc::invalid_geometry);

            for (std::size_t c = 0; c < 6; ++c) {
                const auto part = static_cast<component>(c);
                const auto offset = component_offset(part);
                const auto i_range = axis_range(bounds->min[0], bounds->max[0], nx,
                                                domain.origin[0], domain.spacing[0], offset[0]);
                const auto j_range = axis_range(bounds->min[1], bounds->max[1], ny,
                                                domain.origin[1], domain.spacing[1], offset[1]);
                const auto k_range = axis_range(bounds->min[2], bounds->max[2], nz,
                                                domain.origin[2], domain.spacing[2], offset[2]);
                if (i_range.empty() || j_range.empty() || k_range.empty()) continue;

                auto& ids = map.ids[c];
                for (std::size_t k = k_range.begin; k < k_range.end; ++k) {
                    const double z = domain.origin[2] + (static_cast<double>(k) + offset[2]) * dz;
                    const auto section = common::section_at(solid.shape, z);
                    if (!section) {
                        if (section.error() == common::errc::out_of_range) continue;
                        return fail(errc::invalid_geometry);
                    }
                    for (std::size_t j = j_range.begin; j < j_range.end; ++j) {
                        const double y = domain.origin[1] + (static_cast<double>(j) + offset[1]) * domain.spacing[1];
                        const std::size_t row = k * nx * ny + j * nx;
                        for (std::size_t i = i_range.begin; i < i_range.end; ++i) {
                            const double x =
                                domain.origin[0] + (static_cast<double>(i) + offset[0]) * domain.spacing[0];
                            if (boost::geometry::covered_by(common::point2d{x, y}, *section))
                                ids[row + i] = static_cast<std::uint32_t>(solid.material);
                        }
                    }
                }
            }
        }
        return map;
    });
}

} // namespace litho::fdtd
