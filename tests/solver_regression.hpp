#pragma once

#include <litho_fdtd/fdtd.hpp>

// Shared scenarios for serial/parallel and pre-optimization numerical checks.
// Odd row sizes exercise SIMD tails; transverse planes include CPML cells.
inline litho::fdtd::problem kernel_regression_problem(std::size_t scenario) {
    using namespace litho::fdtd;
    problem p;
    p.grid.cells = scenario == 8 ? std::array<std::size_t, 3>{1, 2, 96}
                                 : std::array<std::size_t, 3>{17, 13, 96};
    p.grid.spacing = {10e-9, 10e-9, 10e-9};
    (void)p.materials.add(material{"air", 1, 0, 1, std::nullopt});
    (void)p.materials.add(material{"lossy magnetic dielectric", 2.25, 1400, 1.3, std::nullopt});
    if (scenario == 6 || scenario == 7)
        (void)p.materials.add(material{"Drude", 1, 0, 1, drude_pole{1.5e16, 2e15}});
    p.boundaries.axes = {boundary_kind::periodic, boundary_kind::periodic, boundary_kind::cpml};
    p.boundaries.cpml.layers = 3;
    if (scenario == 4 || scenario == 5)
        p.boundaries.axes = {boundary_kind::cpml, boundary_kind::cpml, boundary_kind::cpml};
    if (scenario == 7) p.boundaries.axes[0] = boundary_kind::cpml;
    p.source.z = 200e-9;
    p.source.wavelength = 100e-9;
    p.source.pulse_fwhm = 0.3;
    p.source.polarization = 0.6;
    if (scenario == 3 || scenario == 4 || scenario == 7) {
        p.source.theta = 0.25;
        p.source.phi = 0.4;
    }
    if (scenario == 1)
        p.geometry.layers.push_back(planar_layer{440e-9, 620e-9, 1});
    if ((scenario >= 2 && scenario <= 3) || (scenario >= 5 && scenario <= 7)) {
        auto rectangle = common::rectangle(common::box<double, 2>{{30e-9, 20e-9}, {120e-9, 90e-9}});
        auto solid = common::make_extrusion(*rectangle, 440e-9, 620e-9, 1.0);
        p.geometry.solids.push_back(patterned_solid{*solid, scenario >= 6 ? 2u : 1u});
    }
    p.run.steps = 120;
    p.run.workers = 1;
    p.run.chunk_rows = 7;
    for (auto normal : {axis::x, axis::y, axis::z}) {
        const auto size = p.grid.cells[static_cast<std::size_t>(normal)];
        const auto padding = is_periodic(p.boundaries, static_cast<std::size_t>(normal))
                                 ? 0 : p.boundaries.cpml.layers;
        for (auto index : {padding, size / 2, size - padding - 1}) {
            monitor m;
            m.normal = normal;
            m.index = index;
            m.frequencies = {speed_of_light / 100e-9, speed_of_light / 110e-9};
            p.monitors.push_back(m);
        }
    }
    return p;
}
