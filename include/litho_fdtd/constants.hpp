#pragma once

#include <numbers>

namespace litho::fdtd {

// SI constants. The solver uses E in V/m and H in A/m.
inline constexpr double speed_of_light = 299792458.0;
inline constexpr double vacuum_permeability = 4.0e-7 * std::numbers::pi;
inline constexpr double vacuum_permittivity = 1.0 / (vacuum_permeability * speed_of_light * speed_of_light);
inline constexpr double vacuum_impedance = vacuum_permeability * speed_of_light;

} // namespace litho::fdtd
