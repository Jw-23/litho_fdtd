#pragma once

#include <litho_fdtd/boundary.hpp>
#include <litho_fdtd/error.hpp>
#include <litho_fdtd/grid.hpp>
#include <litho_fdtd/monitor.hpp>
#include <litho_fdtd/source.hpp>
#include <litho_fdtd/structure.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace litho::fdtd {

struct run_options {
    double courant = 0.95; // Courant number, must be in (0, 1]
    double duration = 60;  // optical periods when `steps` is unset
    std::optional<std::size_t> steps;
    std::size_t workers = 0; // 0 uses hardware concurrency
};

struct problem {
    grid grid;
    material_library materials;
    structure geometry;
    boundaries boundaries;
    plane_wave_source source;
    std::vector<monitor> monitors;
    run_options run;
};

struct run_result {
    double dt = 0;
    std::size_t steps = 0;
    double elapsed_seconds = 0;
    double frequency = 0;
    std::size_t workers = 1;
    std::array<double, 3> source_wavevector{}; // rad/m
    std::vector<monitor_result> monitors;
};

[[nodiscard]] inline result<void> validate(const run_options& value) noexcept {
    if (!std::isfinite(value.courant) || !(value.courant > 0) || value.courant > 1)
        return fail(errc::invalid_configuration);
    if (!std::isfinite(value.duration) || !(value.duration > 0))
        return fail(errc::invalid_configuration);
    if (value.steps && *value.steps == 0) return fail(errc::invalid_configuration);
    return {};
}

// Run a complete FDTD simulation. The function validates the configuration,
// rasterizes the structure, allocates the field state and executes the time
// loop. All fallible work is reported through `std::expected`.
[[nodiscard]] result<run_result> solve(const problem& value) noexcept;

} // namespace litho::fdtd
