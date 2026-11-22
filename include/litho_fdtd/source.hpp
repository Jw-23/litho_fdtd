#pragma once

#include <litho_fdtd/constants.hpp>
#include <litho_fdtd/error.hpp>
#include <litho_geometry/common.hpp>

#include <cmath>
#include <numbers>

namespace litho::fdtd {

// A plane-wave pulse initialized in a uniform, lossless layer. `z` is the
// initial pulse-center coordinate; the pulse is spatially separated from the
// structure by the caller. `theta` is measured from +z and must be below
// grazing; `from_below` selects +z or -z propagation. `polarization` is 0 for s
// (TE) and pi/2 for p (TM). `pulse_fwhm` is expressed in optical periods at the
// vacuum wavelength.
struct plane_wave_source {
    double z = 0;
    double wavelength = 193e-9;
    double theta = 0;
    double phi = 0;
    double polarization = 0;
    double amplitude = 1; // V/m, peak complex envelope
    double pulse_fwhm = 2;
    bool from_below = true;
    double phase = 0;
};

[[nodiscard]] inline double frequency(const plane_wave_source& value) noexcept {
    return speed_of_light / value.wavelength;
}

[[nodiscard]] inline double period(const plane_wave_source& value) noexcept {
    return value.wavelength / speed_of_light;
}

[[nodiscard]] inline double pulse_sigma(const plane_wave_source& value) noexcept {
    constexpr double fwhm_to_sigma = 0.42466090014400953; // 1 / (2 sqrt(2 ln 2))
    return value.pulse_fwhm * period(value) * fwhm_to_sigma;
}

[[nodiscard]] inline result<void> validate(const plane_wave_source& value) noexcept {
    if (!std::isfinite(value.z)) return fail(errc::invalid_argument);
    if (!std::isfinite(value.wavelength) || !(value.wavelength > 0)) return fail(errc::invalid_argument);
    if (!std::isfinite(value.theta) || value.theta < 0 || value.theta >= std::numbers::pi / 2)
        return fail(errc::invalid_argument);
    if (!std::isfinite(value.phi)) return fail(errc::invalid_argument);
    if (!std::isfinite(value.polarization)) return fail(errc::invalid_argument);
    if (!std::isfinite(value.amplitude) || value.amplitude == 0) return fail(errc::invalid_argument);
    if (!std::isfinite(value.pulse_fwhm) || !(value.pulse_fwhm > 0)) return fail(errc::invalid_argument);
    if (!std::isfinite(value.phase)) return fail(errc::invalid_argument);
    return {};
}

// Unit propagation direction in a medium with refractive index n. The caller
// supplies the index at the source plane; the source must sit in a uniform,
// lossless layer.
[[nodiscard]] inline litho::geometry::common::vector3d propagation_direction(const plane_wave_source& value) noexcept {
    const double st = std::sin(value.theta);
    const double ct = std::cos(value.theta);
    const double sign = value.from_below ? 1.0 : -1.0;
    return {st * std::cos(value.phi), st * std::sin(value.phi), sign * ct};
}

} // namespace litho::fdtd
