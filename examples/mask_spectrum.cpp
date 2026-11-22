#include <litho_fdtd/fdtd.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <string>

namespace {

namespace common = litho::geometry::common;

litho::fdtd::material_library make_materials(double wavelength) {
    litho::fdtd::material_library library;
    (void)library.add(litho::fdtd::material{"air", 1.0, 0.0, 1.0, std::nullopt});
    auto glass = litho::fdtd::from_refractive_index("fused silica", {1.56, 0.0}, wavelength);
    auto absorber = litho::fdtd::from_refractive_index("MoSi absorber", {2.5, 0.6}, wavelength);
    if (!glass || !absorber) return library;
    (void)library.add(*glass);
    (void)library.add(*absorber);
    return library;
}

} // namespace

int main() {
    using namespace litho::fdtd;

    const double wavelength = 193e-9;
    const double spacing = 5e-9;
    const double period_x = 400e-9;
    const double period_y = 20e-9;

    mask_problem problem_;
    problem_.grid.cells = {80, 4, 120};
    problem_.grid.spacing = {spacing, spacing, spacing};
    problem_.materials = make_materials(wavelength);
    if (problem_.materials.size() != 3) {
        std::fprintf(stderr, "failed to build the material library\n");
        return 1;
    }

    problem_.mask.background = 0;                                  // air
    problem_.mask.layers.push_back(planar_layer{0.0, 300e-9, 1});  // silica substrate
    problem_.boundaries.cpml.layers = 20;
    problem_.source.z = 120e-9; // initial pulse center inside the substrate
    problem_.source.wavelength = wavelength;
    problem_.source.theta = 0;
    problem_.source.polarization = 0; // y-polarized (TE against the lines)
    problem_.source.amplitude = 1.0;
    problem_.source.pulse_fwhm = 0.5;
    problem_.run.duration = 20;

    problem_.monitor_z = 450e-9; // in air, above the absorber
    problem_.spectrum.min_order_x = -6;
    problem_.spectrum.max_order_x = 6;
    problem_.spectrum.min_order_y = 0;
    problem_.spectrum.max_order_y = 0;
    problem_.spectrum.background_index = 1.0;

    auto rectangle = common::rectangle(common::box<double, 2>{{0.0, 0.0}, {period_x / 2, period_y}});
    if (!rectangle) {
        std::fprintf(stderr, "failed to build the absorber section\n");
        return 1;
    }
    auto line = common::make_extrusion(*rectangle, 300e-9, 360e-9, 1.0);
    if (!line) {
        std::fprintf(stderr, "failed to build the absorber extrusion\n");
        return 1;
    }
    problem_.mask.solids.push_back(patterned_solid{*line, 2}); // MoSi

    auto result = compute_mask_spectrum(problem_);
    if (!result) {
        std::fprintf(stderr, "mask spectrum failed: %s\n", std::string(message(result.error())).c_str());
        return 1;
    }

    const auto& orders = result->normalized ? *result->normalized : result->spectrum;
    std::printf("lambda = %.1f nm, incident power = %.6e W/m^2, transmitted efficiency = %.4f\n",
                wavelength * 1e9, result->incident_power, result->total_efficiency);
    std::printf("  m       kx (rad/m)      |Ey|       phase (rad)   efficiency\n");
    for (const auto& order : result->spectrum.orders) {
        const auto* normalized = orders.find(order.m, order.n);
        const double magnitude = normalized ? std::abs(normalized->field[1]) : 0.0;
        const double phase = normalized ? std::arg(normalized->field[1]) : 0.0;
        std::printf("  %+2d   %12.5e   %9.5f   %12.5f   %10.5f\n", order.m, order.kx, magnitude, phase,
                    order.efficiency.value_or(0.0));
    }

    std::ofstream csv("litho_fdtd_mask_spectrum.csv");
    if (csv) {
        csv << "m,n,kx,ky,|Ey|,phase,efficiency\n";
        for (const auto& order : result->spectrum.orders) {
            const auto* normalized = orders.find(order.m, order.n);
            const double magnitude = normalized ? std::abs(normalized->field[1]) : 0.0;
            const double phase = normalized ? std::arg(normalized->field[1]) : 0.0;
            csv << order.m << ',' << order.n << ',' << order.kx << ',' << order.ky << ',' << magnitude << ','
                << phase << ',' << order.efficiency.value_or(0.0) << '\n';
        }
        std::printf("wrote litho_fdtd_mask_spectrum.csv\n");
    }
    return 0;
}
