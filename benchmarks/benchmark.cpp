#include <litho_fdtd/fdtd.hpp>

#include <cstdio>
#include <cstdlib>
#include <string>


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
    return 0;
}
