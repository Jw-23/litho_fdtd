#include <litho_fdtd/fdtd.hpp>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <memory>
#include <chrono>


int main(int argc, char** argv) {
    using namespace litho::fdtd;

    problem value;
    value.grid.cells = {120, 120, 240};
    value.grid.spacing = {5e-9, 5e-9, 5e-9};
    (void)value.materials.add(material{"air", 1.0, 0.0, 1.0, std::nullopt});
    value.geometry.background = 0;
    value.boundaries.axes = {boundary_kind::periodic, boundary_kind::periodic, boundary_kind::cpml};
    value.boundaries.cpml.layers = 20;
    value.source.z = 300e-9;
    value.source.wavelength = 193e-9;
    value.source.pulse_fwhm = 1.0;
    value.run.steps = 300;
    value.run.workers = argc > 1 ? static_cast<std::size_t>(std::strtoull(argv[1], nullptr, 10)) : 0;

    // Optional monitor workload and a fresh tracing directory.
    if (argc > 2 && std::string(argv[2]) == "monitors") {
        for (std::size_t index : {100u, 140u, 180u}) {
            monitor detector;
            detector.index = index;
            const double f = speed_of_light / value.source.wavelength;
            detector.frequencies = {0.9 * f, f, 1.1 * f, 1.2 * f};
            value.monitors.push_back(std::move(detector));
        }
    }
    std::unique_ptr<litho::Session> trace;
    if (argc > 3) {
        litho::Options options;
        options.directory = argv[3];
        options.build_id = "fdtd-benchmark-release";
        trace = std::make_unique<litho::Session>(options);
        value.run.trace.session = trace.get();
        value.run.trace.context = trace->stream("benchmark", {}, LT_SOURCE);
        value.run.trace.step_interval = 50;
    }
    const auto wall_start = std::chrono::steady_clock::now();
    auto result = solve(value);
    if (!result) {
        std::fprintf(stderr, "benchmark failed: %s\n", std::string(message(result.error())).c_str());
        return 1;
    }
    const double cells = static_cast<double>(value.grid.size());
    const double cell_updates = cells * static_cast<double>(result->steps);
    std::printf("grid %zu x %zu x %zu, %zu steps, %zu worker(s), %.3f s, %.1f Mcell-updates/s\n",
                value.grid.cells[0], value.grid.cells[1], value.grid.cells[2], result->steps, result->workers,
                result->elapsed_seconds, cell_updates / result->elapsed_seconds / 1e6);
    std::printf("total solve %.3f s\n", std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start).count());
    if (trace) {
        trace->close();
        const auto stats = trace->stats();
        std::printf("trace: %llu events written, %llu dropped, error: %s\n",
            static_cast<unsigned long long>(stats.written_events),
            static_cast<unsigned long long>(stats.dropped_events), stats.error.c_str());
        if (!stats.error.empty() || stats.dropped_events != 0) return 1;
    }
    return 0;
}
