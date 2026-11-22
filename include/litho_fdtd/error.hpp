#pragma once

#include <expected>
#include <new>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace litho::fdtd {

enum class errc {
    invalid_argument = 1,
    invalid_configuration,
    unsupported,
    out_of_range,
    invalid_geometry,
    missing_material,
    allocation_failure,
    numerical_failure,
    unstable,
    source_inside_pml,
    monitor_inside_pml,
    source_plane_not_uniform,
    solve_failed
};

constexpr std::string_view message(errc error) noexcept {
    switch (error) {
    case errc::invalid_argument: return "invalid argument";
    case errc::invalid_configuration: return "invalid solver configuration";
    case errc::unsupported: return "unsupported configuration";
    case errc::out_of_range: return "parameter outside its domain";
    case errc::invalid_geometry: return "invalid geometry";
    case errc::missing_material: return "missing material";
    case errc::allocation_failure: return "memory allocation failed";
    case errc::numerical_failure: return "numerical failure";
    case errc::unstable: return "time step violates the CFL condition";
    case errc::source_inside_pml: return "source plane lies inside a CPML region";
    case errc::monitor_inside_pml: return "monitor plane lies inside a CPML region";
    case errc::source_plane_not_uniform: return "source plane crosses a material boundary";
    case errc::solve_failed: return "solver did not produce a result";
    }
    return "unknown FDTD error";
}

constexpr std::unexpected<errc> fail(errc error) noexcept { return std::unexpected(error); }

template<class T>
using result = std::expected<T, errc>;

namespace detail {
// Allocating public entry points translate bad_alloc/length_error once.
template<class F>
[[nodiscard]] auto guard(F&& function) noexcept -> std::invoke_result_t<F> {
    try {
        return std::invoke(std::forward<F>(function));
    } catch (const std::bad_alloc&) {
        return fail(errc::allocation_failure);
    } catch (const std::length_error&) {
        return fail(errc::out_of_range);
    }
}
} // namespace detail

} // namespace litho::fdtd
