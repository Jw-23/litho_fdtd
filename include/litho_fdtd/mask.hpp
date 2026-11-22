#pragma once

#include <litho_fdtd/error.hpp>
#include <litho_fdtd/solver.hpp>
#include <litho_fdtd/spectrum.hpp>

#include <array>
#include <optional>

namespace litho::fdtd {

// High-level mask spectrum request. The monitor is a z-normal plane; the
// periodic x/y directions are the mask period. Set `spectrum.background_index`
// to the refractive index at the monitor plane (1 for air, 1.56 for fused
// silica at 193 nm, and so on).
struct mask_problem {
    grid grid;
    material_library materials;
    structure mask;
    plane_wave_source source;
    boundaries boundaries{{boundary_kind::periodic, boundary_kind::periodic, boundary_kind::cpml}, {}};
    double monitor_z = 0;
    std::optional<std::array<double, 2>> monitor_time_window;
    spectrum_options spectrum;
    run_options run;
    bool use_reference = true;
};

struct mask_spectrum {
    spectrum_result spectrum;
    std::optional<spectrum_result> normalized; // incident-normalized complex orders
    std::optional<spectrum_result> reference;
    double incident_power = 0;
    double total_efficiency = 0;
    double dt = 0;
    std::size_t steps = 0;
    std::array<double, 3> source_wavevector{};
};

// Run the patterned mask and, unless disabled, an unpatterned reference stack
// to obtain the incident power. The returned orders carry efficiencies relative
// to that reference. Both runs use the same source and monitor time window.
[[nodiscard]] result<mask_spectrum> compute_mask_spectrum(const mask_problem& value) noexcept;

} // namespace litho::fdtd
