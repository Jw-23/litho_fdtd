#pragma once

#include <litho_fdtd/constants.hpp>
#include <litho_fdtd/error.hpp>

#include <cmath>
#include <complex>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace litho::fdtd {

// Time convention: fields vary as exp(-i omega t). A complex refractive index
// is therefore n = n' + i n'', with n'' >= 0 for a passive material.
struct drude_pole {
    double plasma_frequency = 0; // rad/s
    double collision_rate = 0;   // rad/s
};

// Isotropic material. Non-dispersive loss is represented by a conductivity;
// a single Drude pole is used when the real part of the permittivity at the
// reference wavelength is negative, which keeps the time-domain update stable
// for metals. Both models reproduce the requested complex permittivity exactly
// at the reference angular frequency.
struct material {
    std::string name;
    double epsilon_inf = 1;
    double conductivity = 0; // S/m
    double mu_r = 1;
    std::optional<drude_pole> drude;

    [[nodiscard]] bool dispersive() const noexcept { return drude.has_value(); }

    [[nodiscard]] std::complex<double> relative_permittivity(double omega) const noexcept {
        std::complex<double> value = epsilon_inf;
        if (conductivity != 0 && omega != 0)
            value += std::complex<double>(0.0, conductivity / (omega * vacuum_permittivity));
        if (drude) {
            const double wp2 = drude->plasma_frequency * drude->plasma_frequency;
            value -= wp2 / std::complex<double>(omega * omega, drude->collision_rate * omega);
        }
        return value;
    }
};

[[nodiscard]] inline result<void> validate(const material& value) noexcept {
    if (!std::isfinite(value.epsilon_inf) || !(value.epsilon_inf > 0))
        return fail(errc::invalid_argument);
    if (!std::isfinite(value.conductivity) || value.conductivity < 0)
        return fail(errc::invalid_argument);
    if (!std::isfinite(value.mu_r) || !(value.mu_r > 0))
        return fail(errc::invalid_argument);
    if (value.drude) {
        if (!std::isfinite(value.drude->plasma_frequency) || !(value.drude->plasma_frequency > 0))
            return fail(errc::invalid_argument);
        if (!std::isfinite(value.drude->collision_rate) || !(value.drude->collision_rate > 0))
            return fail(errc::unsupported);
    }
    return {};
}

namespace detail {
[[nodiscard]] inline result<material> from_complex_permittivity(std::string name,
                                                               std::complex<double> epsilon,
                                                               double wavelength) noexcept {
    if (!std::isfinite(wavelength) || !(wavelength > 0)) return fail(errc::invalid_argument);
    if (!std::isfinite(epsilon.real()) || !std::isfinite(epsilon.imag()))
        return fail(errc::invalid_argument);
    if (epsilon.imag() < 0) return fail(errc::invalid_argument); // active material
    if (epsilon.real() == 0) return fail(errc::invalid_argument);

    const double omega = 2.0 * std::numbers::pi * speed_of_light / wavelength;
    if (epsilon.real() > 0) {
        material value{std::move(name), epsilon.real(), epsilon.imag() * omega * vacuum_permittivity, 1, std::nullopt};
        if (auto valid = validate(value); !valid) return fail(valid.error());
        return value;
    }

    // Drude fit at omega. With epsilon_inf = 1:
    //   epsilon - 1 = -wp^2 / (omega^2 + i gamma omega).
    const std::complex<double> delta = epsilon - std::complex<double>(1.0, 0.0);
    const double a = -delta.real();
    const double b = delta.imag();
    if (!(a > 0) || !(b > 0)) return fail(errc::unsupported);
    const double ratio = b / a;
    const double gamma = omega * ratio;
    const double wp2 = a * (omega * omega + gamma * gamma);
    if (!std::isfinite(gamma) || !(gamma > 0) || !std::isfinite(wp2) || !(wp2 > 0))
        return fail(errc::numerical_failure);

    material value{std::move(name), 1.0, 0.0, 1.0, drude_pole{std::sqrt(wp2), gamma}};
    if (auto valid = validate(value); !valid) return fail(valid.error());
    return value;
}
} // namespace detail

[[nodiscard]] inline result<material> from_refractive_index(std::string name,
                                                            std::complex<double> n,
                                                            double wavelength) noexcept {
    if (n.imag() < 0) return fail(errc::invalid_argument);
    return detail::from_complex_permittivity(std::move(name), n * n, wavelength);
}

[[nodiscard]] inline result<material> from_permittivity(std::string name,
                                                        std::complex<double> epsilon,
                                                        double wavelength) noexcept {
    return detail::from_complex_permittivity(std::move(name), epsilon, wavelength);
}

struct material_library {
    std::vector<material> entries;

    [[nodiscard]] result<std::size_t> add(material value) noexcept {
        return detail::guard([&]() -> result<std::size_t> {
            if (auto valid = validate(value); !valid) return fail(valid.error());
            entries.push_back(std::move(value));
            return entries.size() - 1;
        });
    }
    [[nodiscard]] const material* at(std::size_t index) const noexcept {
        return index < entries.size() ? &entries[index] : nullptr;
    }
    [[nodiscard]] std::span<const material> all() const noexcept { return entries; }
    [[nodiscard]] std::size_t size() const noexcept { return entries.size(); }
};

} // namespace litho::fdtd
