#include "test_support.hpp"

#include <litho_fdtd/fdtd.hpp>

#include <cmath>
#include <complex>
#include <numbers>

using namespace litho::fdtd;
namespace common = litho::geometry::common;

namespace {

constexpr double wavelength = 200e-9;
constexpr double k0 = 2.0 * std::numbers::pi / wavelength;

problem vacuum_problem(std::size_t nx, std::size_t ny, std::size_t nz, double spacing) {
    problem value;
    value.grid.cells = {nx, ny, nz};
    value.grid.spacing = {spacing, spacing, spacing};
    CHECK(value.materials.add(material{"air", 1.0, 0.0, 1.0, std::nullopt}).has_value());
    value.geometry.background = 0;
    value.boundaries.axes = {boundary_kind::periodic, boundary_kind::periodic, boundary_kind::cpml};
    value.boundaries.cpml.layers = 20;
    value.source.wavelength = wavelength;
    value.source.amplitude = 1.0;
    return value;
}

void plane_wave_vacuum() {
    auto value = vacuum_problem(4, 4, 120, 10e-9);
    value.source.z = 250e-9;
    value.source.polarization = 0;
    value.source.pulse_fwhm = 1.0;
    value.run.duration = 10;
    value.run.workers = 4;

    monitor detector;
    detector.name = "vacuum";
    detector.normal = axis::z;
    detector.index = 60;
    spectrum_options options;
    options.min_order_x = options.max_order_x = 0;
    options.min_order_y = options.max_order_y = 0;
    options.background_index = 1.0;
    detector.spectrum = options;
    value.monitors.push_back(detector);

    auto run = solve(value);
    CHECK(run.has_value());
    if (!run) return;
    CHECK(run->workers == 4);
    const auto& monitor = run->monitors.at(0);
    CHECK(monitor.spectrum.has_value());
    const auto* order = monitor.spectrum ? monitor.spectrum->find(0, 0) : nullptr;
    CHECK(order != nullptr);
    if (!order) return;

    CHECK_NEAR(std::abs(order->field[1]), 1.0, 0.06);
    CHECK_NEAR(order->power, 0.5 / vacuum_impedance, 0.06);

    // For a +z wave with E along y the magnetic field is -x.
    const auto hx = monitor.component_values(0, component::hx);
    std::complex<double> average{};
    for (const auto field : hx) average += field;
    average /= static_cast<double>(hx.size());
    const double ratio = (average / order->field[1]).real();
    CHECK_NEAR(ratio, -1.0 / vacuum_impedance, 0.03);
}

void fresnel_reflection() {
    const double spacing = 5e-9;
    auto base = vacuum_problem(4, 4, 280, spacing);
    CHECK(base.materials.add(material{"glass", 4.0, 0.0, 1.0, std::nullopt}).has_value());
    base.boundaries.cpml.layers = 40;
    base.source.z = 250e-9;
    base.source.polarization = 0;
    base.source.pulse_fwhm = 0.5;
    base.run.duration = 10;

    monitor detector;
    detector.name = "reflection";
    detector.normal = axis::z;
    detector.index = 120;
    spectrum_options options;
    options.min_order_x = options.max_order_x = 0;
    options.min_order_y = options.max_order_y = 0;
    options.background_index = 1.0;
    detector.spectrum = options;
    const double period = wavelength / speed_of_light;
    detector.time_window = std::array<double, 2>{0.0, 5.0 * period};
    base.monitors.push_back(detector);

    auto reference = base;
    auto full = base;
    const double interface_z = 700e-9;
    full.geometry.layers.push_back(planar_layer{interface_z, 1400e-9, 1});

    auto reference_run = solve(reference);
    auto full_run = solve(full);
    CHECK(reference_run.has_value());
    CHECK(full_run.has_value());
    if (!reference_run || !full_run) return;

    const auto* incident = reference_run->monitors[0].spectrum->find(0, 0);
    const auto* total = full_run->monitors[0].spectrum->find(0, 0);
    CHECK(incident && total);
    if (!incident || !total) return;
    const std::complex<double> ratio = (total->field[1] - incident->field[1]) / incident->field[1];
    const std::complex<double> measured =
        ratio * std::exp(std::complex<double>(0.0, -2.0 * k0 * (interface_z - 600e-9)));
    const std::complex<double> expected((1.0 - 2.0) / (1.0 + 2.0), 0.0);
    CHECK(std::abs(measured - expected) < 0.08);
}

void oblique_bloch(std::complex<double> polarization, bool check_transverse) {
    auto value = vacuum_problem(80, 8, 240, 5e-9);
    value.boundaries.cpml.layers = 40;
    value.source.z = 250e-9;
    value.source.theta = 20.0 * std::numbers::pi / 180.0;
    value.source.polarization = polarization.real();
    value.source.pulse_fwhm = 1.0;
    value.run.duration = 10;

    monitor detector;
    detector.name = "oblique";
    detector.normal = axis::z;
    detector.index = 90;
    spectrum_options options;
    options.min_order_x = -1;
    options.max_order_x = 1;
    options.min_order_y = options.max_order_y = 0;
    options.background_index = 1.0;
    detector.spectrum = options;
    value.monitors.push_back(detector);

    auto run = solve(value);
    CHECK(run.has_value());
    if (!run) return;
    CHECK_NEAR(run->source_wavevector[0], k0 * std::sin(value.source.theta), 1e-6);
    CHECK_NEAR(run->source_wavevector[2], k0 * std::cos(value.source.theta), 1e-6);

    const auto* order = run->monitors[0].spectrum ? run->monitors[0].spectrum->find(0, 0) : nullptr;
    CHECK(order != nullptr);
    if (!order) return;
    const double amplitude = std::sqrt(std::norm(order->field[0]) + std::norm(order->field[1]) +
                                       std::norm(order->field[2]));
    CHECK_NEAR(amplitude, 1.0, 0.06);
    CHECK_NEAR(order->power, 0.5 / vacuum_impedance * std::cos(value.source.theta), 0.08);

    if (check_transverse) {
        const std::complex<double> dot =
            order->field[0] * order->kx + order->field[1] * order->ky + order->field[2] * order->kz;
        CHECK(std::abs(dot) < 0.12 * amplitude * std::abs(order->kz));
    }

    if (run->monitors[0].spectrum) {
        for (const auto& candidate : run->monitors[0].spectrum->orders)
            if (candidate.m != 0 || candidate.n != 0) CHECK(candidate.power < 1e-6);
    }
}

monitor_result subtract(const monitor_result& full, const monitor_result& reference) {
    monitor_result out = full;
    for (std::size_t i = 0; i < out.values.size(); ++i) out.values[i] -= reference.values[i];
    out.spectrum.reset();
    return out;
}

void grating_energy() {
    problem base;
    base.grid.cells = {30, 4, 120};
    base.grid.spacing = {10e-9, 10e-9, 10e-9};
    CHECK(base.materials.add(material{"air", 1.0, 0.0, 1.0, std::nullopt}).has_value());
    CHECK(base.materials.add(material{"glass", 4.0, 0.0, 1.0, std::nullopt}).has_value());
    base.geometry.background = 0;
    base.boundaries.axes = {boundary_kind::periodic, boundary_kind::periodic, boundary_kind::cpml};
    base.boundaries.cpml.layers = 20;
    base.source.z = 250e-9;
    base.source.wavelength = wavelength;
    base.source.polarization = 0;
    base.source.pulse_fwhm = 0.5;
    base.run.duration = 20;

    spectrum_options options;
    options.min_order_x = -2;
    options.max_order_x = 2;
    options.min_order_y = options.max_order_y = 0;
    options.background_index = 1.0;

    monitor transmitted;
    transmitted.name = "transmitted";
    transmitted.normal = axis::z;
    transmitted.index = 90;
    transmitted.spectrum = options;
    monitor reflected;
    reflected.name = "reflected";
    reflected.normal = axis::z;
    reflected.index = 60;
    reflected.spectrum = options;
    base.monitors = {transmitted, reflected};

    auto rectangle = common::rectangle(common::box<double, 2>{{0, 0}, {150e-9, 40e-9}});
    CHECK(rectangle.has_value());
    auto solid = common::make_extrusion(*rectangle, 700e-9, 800e-9, 1.0);
    CHECK(solid.has_value());
    auto full = base;
    full.geometry.solids.push_back(patterned_solid{*solid, 1});

    auto reference_run = solve(base);
    CHECK(reference_run.has_value());
    if (!reference_run) return;
    const auto* incident = reference_run->monitors[0].spectrum->find(0, 0);
    CHECK(incident != nullptr);
    if (!incident) return;
    const double incident_power = incident->power;
    CHECK(incident_power > 0);

    full.monitors[0].spectrum->incident_power = incident_power;
    auto full_run = solve(full);
    CHECK(full_run.has_value());
    if (!full_run) return;
    const double transmitted_efficiency = full_run->monitors[0].spectrum->total_efficiency.value_or(0);

    auto reflected_difference = subtract(full_run->monitors[1], reference_run->monitors[1]);
    auto reflected_options = options;
    reflected_options.incident_power = incident_power;
    reflected_options.propagation_sign = -1;
    auto reflected_spectrum = compute_spectrum(reflected_difference, reflected_options);
    CHECK(reflected_spectrum.has_value());
    if (!reflected_spectrum) return;
    const double reflected_efficiency = reflected_spectrum->total_efficiency.value_or(0);

    CHECK_NEAR(transmitted_efficiency, 0.54, 0.06);
    CHECK_NEAR(reflected_efficiency, 0.45, 0.06);
    CHECK_NEAR(transmitted_efficiency + reflected_efficiency, 1.0, 0.03);
}

void mask_spectrum_api() {
    mask_problem value;
    value.grid.cells = {30, 4, 120};
    value.grid.spacing = {10e-9, 10e-9, 10e-9};
    CHECK(value.materials.add(material{"air", 1.0, 0.0, 1.0, std::nullopt}).has_value());
    CHECK(value.materials.add(material{"glass", 4.0, 0.0, 1.0, std::nullopt}).has_value());
    value.mask.background = 0;
    value.boundaries.cpml.layers = 20;
    value.source.z = 250e-9;
    value.source.wavelength = wavelength;
    value.source.polarization = 0;
    value.source.pulse_fwhm = 0.5;
    value.run.duration = 20;
    value.monitor_z = 900e-9;
    value.spectrum.min_order_x = -2;
    value.spectrum.max_order_x = 2;
    value.spectrum.min_order_y = value.spectrum.max_order_y = 0;
    value.spectrum.background_index = 1.0;

    auto rectangle = common::rectangle(common::box<double, 2>{{0, 0}, {150e-9, 40e-9}});
    CHECK(rectangle.has_value());
    auto solid = common::make_extrusion(*rectangle, 700e-9, 800e-9, 1.0);
    CHECK(solid.has_value());
    value.mask.solids.push_back(patterned_solid{*solid, 1});

    auto result = compute_mask_spectrum(value);
    CHECK(result.has_value());
    if (!result) return;
    CHECK(result->incident_power > 0);
    CHECK_NEAR(result->total_efficiency, 0.54, 0.06);
    CHECK(result->normalized.has_value());
    CHECK_NEAR(result->spectrum.total_efficiency.value_or(0), result->total_efficiency, 1e-12);
}

void drude_absorber() {
    auto metal = from_refractive_index("metal", {0.5, 2.0}, wavelength);
    CHECK(metal.has_value());
    if (!metal) return;

    problem value;
    value.grid.cells = {8, 4, 280};
    value.grid.spacing = {5e-9, 5e-9, 5e-9};
    CHECK(value.materials.add(material{"air", 1.0, 0.0, 1.0, std::nullopt}).has_value());
    CHECK(value.materials.add(*metal).has_value());
    value.geometry.background = 0;
    value.geometry.layers.push_back(planar_layer{800e-9, 850e-9, 1});
    value.boundaries.axes = {boundary_kind::periodic, boundary_kind::periodic, boundary_kind::cpml};
    value.boundaries.cpml.layers = 40;
    value.source.z = 250e-9;
    value.source.wavelength = wavelength;
    value.source.polarization = 0;
    value.source.pulse_fwhm = 0.5;
    value.run.duration = 15;

    monitor detector;
    detector.name = "transmitted";
    detector.normal = axis::z;
    detector.index = 180;
    spectrum_options options;
    options.min_order_x = options.max_order_x = 0;
    options.min_order_y = options.max_order_y = 0;
    options.background_index = 1.0;
    detector.spectrum = options;
    value.monitors.push_back(detector);

    auto run = solve(value);
    CHECK(run.has_value());
    if (!run) return;
    const auto* order = run->monitors[0].spectrum->find(0, 0);
    CHECK(order != nullptr);
    if (!order) return;
    CHECK(std::isfinite(order->field[1].real()) && std::isfinite(order->field[1].imag()));
    CHECK(order->power < 1e-4);
}

} // namespace

void run_solver_tests() {
    plane_wave_vacuum();
    fresnel_reflection();
    oblique_bloch(0.0, false);
    oblique_bloch(std::numbers::pi / 2.0, true);
    grating_energy();
    mask_spectrum_api();
    drude_absorber();
}
