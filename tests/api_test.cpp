#include "test_support.hpp"
#include <litho_fdtd/fdtd.hpp>
#include <chrono>
#include <filesystem>
#include <limits>
#include <map>
#include <numbers>

using namespace litho::fdtd;
namespace {
void api_checks() {
    grid domain;
    domain.cells = {4, 5, 6};
    domain.spacing = {0.1, 0.2, 0.3};
    domain.origin = {-0.4, 0.5, -0.6};
    auto plane = monitor_at(domain, axis::z, 0.3);
    CHECK(plane && plane->index == 3 && plane->normal == axis::z);
    CHECK(!monitor_at(domain, axis::z, 0.31));
    CHECK(!monitor_at(domain, axis::z, 1e300));
    CHECK(!monitor_at(domain, static_cast<axis>(5), 0));
    monitor invalid;
    invalid.normal = static_cast<axis>(3);
    CHECK(!validate(invalid, domain));
    monitor implicit_frequency;
    implicit_frequency.spectrum = spectrum_options{};
    implicit_frequency.spectrum->frequency_index = 1;
    CHECK(!validate(implicit_frequency, domain));
    monitor_result truncated;
    truncated.shape = {4, 5};
    truncated.frequencies = {1};
    truncated.values.resize(3);
    CHECK(truncated.component_values(0, component::ex).empty());
    truncated.values.resize(120);
    CHECK(truncated.component_values(0, static_cast<component>(6)).empty());
    CHECK(truncated.component_values(std::numeric_limits<std::size_t>::max(), component::ex).empty());
    run_options options;
    options.steps = 1;
    options.duration = std::numeric_limits<double>::quiet_NaN();
    CHECK(validate(options).has_value());
    options.steps.reset();
    CHECK(!validate(options));
}

void spectrum_regression() {
    // Independent direct DFT using the pre-optimization sampling convention.
    monitor_result plane;
    plane.shape = {7, 5};
    plane.spacing = {7e-9, 9e-9};
    plane.origin = {-13e-9, 21e-9};
    plane.normal_spacing = 8e-9;
    plane.frequencies = {speed_of_light / 200e-9, speed_of_light / 193e-9};
    plane.values.resize(12 * plane.cell_count());
    for (std::size_t i = 0; i < plane.values.size(); ++i)
        plane.values[i] = {std::sin(0.7 * i), std::cos(0.13 * i)};
    spectrum_options options;
    options.frequency_index = 1;
    options.min_order_x = options.min_order_y = -2;
    options.max_order_x = options.max_order_y = 2;
    options.incident_wavevector = std::array<double, 2>{1.2e6, -2.3e6};
    auto spectrum = compute_spectrum(plane, options);
    CHECK(spectrum.has_value());
    if (!spectrum) return;
    for (const auto& order : spectrum->orders) {
        std::array<std::complex<double>, 5> sums{};
        for (std::size_t c = 0; c < 5; ++c) {
            const bool offset_x = c == 0 || c == 2 || c == 4;
            const bool offset_y = c == 1 || c == 2 || c == 3;
            for (std::size_t j = 0; j < plane.shape[1]; ++j)
                for (std::size_t i = 0; i < plane.shape[0]; ++i) {
                    const double x = plane.origin[0] + (i + 0.5 * offset_x) * plane.spacing[0];
                    const double y = plane.origin[1] + (j + 0.5 * offset_y) * plane.spacing[1];
                    sums[c] += plane.at(1, static_cast<component>(c), i + 7 * j) *
                        std::exp(std::complex<double>{0, -order.kx * x - order.ky * y});
                }
            sums[c] /= plane.cell_count();
        }
        for (std::size_t c = 0; c < 3; ++c)
            CHECK(std::abs(order.field[c] - sums[c]) < 2e-14);
        if (order.propagating) {
            const auto stagger = std::exp(std::complex<double>{0, -0.5 * plane.normal_spacing} * order.kz);
            const double power = 0.5 * (sums[0] * std::conj(sums[4] * stagger) -
                                         sums[1] * std::conj(sums[3] * stagger)).real();
            CHECK(std::abs(order.power - power) < 2e-14);
        }
    }
    plane.values.clear();
    CHECK(!compute_spectrum(plane, options));
}

void trace_checks() {
    problem p;
    p.grid.cells = {4, 4, 120};
    p.grid.spacing = {10e-9, 10e-9, 10e-9};
    CHECK(p.materials.add(material{"air", 1, 0, 1, std::nullopt}).has_value());
    p.boundaries.axes = {boundary_kind::periodic, boundary_kind::periodic, boundary_kind::cpml};
    p.boundaries.cpml.layers = 20;
    p.source.z = 300e-9;
    p.source.wavelength = 200e-9;
    p.source.pulse_fwhm = 0.5;
    p.run.steps = 5;
    p.run.workers = 1;
    for (auto normal : {axis::x, axis::y, axis::z}) {
        monitor m;
        m.normal = normal;
        m.index = normal == axis::z ? 60 : 1;
        m.frequencies = {speed_of_light / 200e-9, speed_of_light / 210e-9};
        p.monitors.push_back(m);
    }
    CHECK(validate(p).has_value());
    const auto plain = solve(p);
    CHECK(plain.has_value());
    const auto directory = std::filesystem::temp_directory_path() / ("fdtd-trace-test-" +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    if (plain) {
        // In this transverse-uniform vacuum wave, x/y planes agree, and their
        // z=60 row equals every sample of the z-normal monitor.
        CHECK(plain->monitors[0].values == plain->monitors[1].values);
        for (std::size_t f = 0; f < 2; ++f)
            for (std::size_t c = 0; c < 6; ++c) {
                const auto component = static_cast<litho::fdtd::component>(c);
                const auto sample = plain->monitors[0].at(f, component, 4 * 60);
                for (const auto z_sample : plain->monitors[2].component_values(f, component))
                    CHECK(std::abs(sample - z_sample) < 1e-14);
            }
    }
    litho::Options options;
    options.directory = directory.string();
    litho::Session session(options);
    const auto root = session.stream("test");
    p.run.trace = {&session, root, 2};
    const auto traced = solve(p);
    CHECK(traced.has_value());
    if (plain && traced) {
        for (std::size_t m = 0; m < 3; ++m)
            CHECK(plain->monitors[m].values == traced->monitors[m].values);
    }
    p.run.trace.step_interval = 0;
    CHECK(solve(p).has_value());
    mask_problem mask;
    mask.grid = p.grid;
    mask.materials = p.materials;
    mask.boundaries = p.boundaries;
    mask.source = p.source;
    mask.monitor_z = 600e-9;
    mask.run = p.run;
    mask.run.steps = 120;
    mask.spectrum.min_order_x = mask.spectrum.max_order_x = 0;
    mask.spectrum.min_order_y = mask.spectrum.max_order_y = 0;
    CHECK(compute_mask_spectrum(mask).has_value());
    p.run.steps = 0;
    CHECK(!solve(p));
    session.close();
    CHECK(session.stats().dropped_events == 0);
    CHECK(session.stats().error.empty());
    const auto events = litho::read_events((directory / "events.lt").string());
    CHECK(events.error.empty() && !events.incomplete_tail);
    std::map<std::uint64_t, std::uint64_t> parents;
    std::size_t begins = 0, ends = 0, progress = 0;
    bool failed = false, prepared = false, reference = false, patterned = false;
    for (const auto& event : events.events) {
        if (event.kind == litho::EventKind::Begin) {
            ++begins;
            parents[event.context.span_id] = event.context.parent_span_id;
            CHECK(event.context.stream_id == root.stream_id);
            CHECK(!event.file.empty() && event.line != 0);
        }
        if (event.kind == litho::EventKind::End) ++ends;
        if (event.name == "progress") ++progress;
        failed |= event.name == "failed";
        prepared |= event.name == "prepared";
        reference |= event.name == "reference";
        patterned |= event.name == "patterned";
    }
    CHECK(begins == ends && begins >= 10);
    CHECK(progress == 3 && failed && prepared && reference && patterned);
    for (const auto& [id, parent] : parents) {
        (void)id;
        CHECK(parent == root.span_id || parents.contains(parent));
    }
    // A closed session must not affect numerical work or terminate noexcept APIs.
    p.run.steps = 1;
    CHECK(solve(p).has_value());
    std::filesystem::remove_all(directory);
}
} // namespace
void run_api_tests() { api_checks(); spectrum_regression(); trace_checks(); }
