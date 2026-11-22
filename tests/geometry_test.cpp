#include "test_support.hpp"

#include <litho_fdtd/structure.hpp>
#include <litho_geometry/measure.hpp>

#include <boost/geometry.hpp>

namespace {

namespace common = litho::geometry::common;
namespace gmeasure = litho::geometry::measure;
using namespace litho::fdtd;

void extrusion_slices_and_bounds() {
    auto rectangle = common::rectangle(common::box<double, 2>{{0, 0}, {10, 10}});
    CHECK(rectangle.has_value());
    auto solid = common::make_extrusion(*rectangle, 0.0, 10.0, 2.0);
    CHECK(solid.has_value());
    CHECK_NEAR(gmeasure::area(*rectangle), 100.0, 1e-12);

    auto bottom = common::section_at(*solid, 0.0);
    auto middle = common::section_at(*solid, 5.0);
    auto top = common::section_at(*solid, 10.0);
    CHECK(bottom.has_value() && middle.has_value() && top.has_value());
    CHECK_NEAR(gmeasure::area(*bottom), 100.0, 1e-12);
    CHECK_NEAR(gmeasure::area(*middle), 225.0, 1e-12);
    CHECK_NEAR(gmeasure::area(*top), 400.0, 1e-12);

    auto outside = common::section_at(*solid, 11.0);
    CHECK(!outside.has_value());
    CHECK(outside.error() == common::errc::out_of_range);

    auto bounds = common::bounds(*solid);
    CHECK(bounds.has_value());
    CHECK_NEAR(bounds->min[0], -5.0, 1e-12);
    CHECK_NEAR(bounds->max[0], 15.0, 1e-12);
    CHECK_NEAR(bounds->min[2], 0.0, 1e-12);
    CHECK_NEAR(bounds->max[2], 10.0, 1e-12);

    // The variant API must preserve the same slice semantics.
    common::solid3d variant{*solid};
    auto variant_section = common::section_at(variant, 10.0);
    CHECK(variant_section.has_value());
    CHECK_NEAR(gmeasure::area(*variant_section), 400.0, 1e-12);
    auto variant_bounds = common::bounds(variant);
    CHECK(variant_bounds.has_value());
    CHECK_NEAR(variant_bounds->max[1], 15.0, 1e-12);

    auto bad_scale = common::make_extrusion(*rectangle, 0.0, 10.0, 0.0);
    CHECK(!bad_scale.has_value());
    auto bad_box = common::make_box3<double>({1, 1, 1}, {1, 2, 2});
    CHECK(!bad_box.has_value());
}

void voxelization() {
    grid domain;
    domain.cells = {10, 10, 10};
    domain.spacing = {1.0, 1.0, 1.0};

    material_library materials;
    CHECK(materials.add(material{"air", 1.0, 0.0, 1.0, std::nullopt}).has_value());
    CHECK(materials.add(material{"glass", 4.0, 0.0, 1.0, std::nullopt}).has_value());

    auto rectangle = common::rectangle(common::box<double, 2>{{2, 2}, {8, 8}});
    CHECK(rectangle.has_value());
    auto solid = common::make_extrusion(*rectangle, 2.0, 8.0, 1.0);
    CHECK(solid.has_value());

    structure layout;
    layout.background = 0;
    layout.solids.push_back(patterned_solid{*solid, 1});
    auto map = voxelize(layout, materials, domain);
    CHECK(map.has_value());

    const auto inside = domain.index(5, 5, 5);
    const auto outside = domain.index(0, 5, 5);
    CHECK(map->at(component::ex, inside) == 1);
    CHECK(map->at(component::ey, inside) == 1);
    CHECK(map->at(component::ez, inside) == 1);
    CHECK(map->at(component::ex, outside) == 0);
    CHECK(map->at(component::hx, inside) == 1);

    structure box_structure;
    box_structure.background = 0;
    box_structure.solids.push_back(patterned_solid{*common::make_box3<double>({2, 2, 2}, {8, 8, 8}), 1});
    auto box_map = voxelize(box_structure, materials, domain);
    CHECK(box_map.has_value());
    CHECK(box_map->at(component::ex, inside) == 1);
    CHECK(box_map->at(component::ex, outside) == 0);

    // A layer fills the whole xy plane, and a later solid overrides it.
    structure layered;
    layered.background = 0;
    layered.layers.push_back(planar_layer{0.0, 10.0, 1});
    auto layered_map = voxelize(layered, materials, domain);
    CHECK(layered_map.has_value());
    CHECK(layered_map->at(component::ex, outside) == 1);
}

} // namespace

void run_geometry_tests() {
    extrusion_slices_and_bounds();
    voxelization();
}
