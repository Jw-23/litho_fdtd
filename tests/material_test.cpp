#include "test_support.hpp"

#include <litho_fdtd/material.hpp>

#include <cmath>
#include <complex>
#include <numbers>

using namespace litho::fdtd;

namespace {

constexpr double wavelength = 200e-9;
const double omega = 2.0 * std::numbers::pi * speed_of_light / wavelength;

void transparent_dielectric() {
    auto glass = from_refractive_index("glass", {1.5, 0.0}, wavelength);
    CHECK(glass.has_value());
    CHECK(!glass->dispersive());
    CHECK_NEAR(glass->epsilon_inf, 2.25, 1e-12);
    CHECK_NEAR(glass->conductivity, 0.0, 1e-12);
    const auto epsilon = glass->relative_permittivity(omega);
    CHECK_NEAR(epsilon.real(), 2.25, 1e-12);
    CHECK_NEAR(epsilon.imag(), 0.0, 1e-12);
}

void lossy_dielectric() {
    // n = 1.5 + 0.1i => eps = 2.24 + 0.3i with the exp(-i omega t) convention.
    auto absorber = from_refractive_index("absorber", {1.5, 0.1}, wavelength);
    CHECK(absorber.has_value());
    CHECK_NEAR(absorber->epsilon_inf, 2.24, 1e-12);
    CHECK_NEAR(absorber->conductivity, omega * vacuum_permittivity * 0.3, 1e-10);
    const auto epsilon = absorber->relative_permittivity(omega);
    CHECK_COMPLEX_NEAR(epsilon, std::complex<double>(2.24, 0.3), 1e-10);
}

void drude_metal() {
    // A negative real permittivity is fitted with a Drude pole, not a negative
    // static epsilon, so the time-domain update stays passive and stable.
    const std::complex<double> target(-3.75, 2.0);
    auto metal = from_permittivity("metal", target, wavelength);
    CHECK(metal.has_value());
    CHECK(metal->dispersive());
    CHECK_NEAR(metal->epsilon_inf, 1.0, 1e-12);
    const auto epsilon = metal->relative_permittivity(omega);
    CHECK_COMPLEX_NEAR(epsilon, target, 1e-10);
}

void rejects_active_and_invalid_materials() {
    auto active = from_refractive_index("active", {1.5, -0.1}, wavelength);
    CHECK(!active.has_value());
    auto zero = from_permittivity("zero", {0.0, 0.0}, wavelength);
    CHECK(!zero.has_value());
    auto negative_mu = material{"bad", 1.0, 0.0, -1.0, std::nullopt};
    CHECK(!validate(negative_mu).has_value());

    material_library library;
    CHECK(library.add(material{"air", 1.0, 0.0, 1.0, std::nullopt}).has_value());
    CHECK(library.size() == 1);
    CHECK(library.at(0) != nullptr);
    CHECK(library.at(1) == nullptr);
}

} // namespace

void run_material_tests() {
    transparent_dielectric();
    lossy_dielectric();
    drude_metal();
    rejects_active_and_invalid_materials();
}
