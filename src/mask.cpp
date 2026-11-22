#include <litho_fdtd/mask.hpp>

#include <cmath>
#include <utility>

namespace litho::fdtd {

namespace {

result<std::size_t> plane_index(const grid& domain, double z) noexcept {
    if (!std::isfinite(z)) return fail(errc::invalid_argument);
    const double position = (z - domain.origin[2]) / domain.spacing[2];
    const auto nearest = static_cast<long long>(std::llround(position));
    if (std::abs(position - static_cast<double>(nearest)) > 1e-6)
        return fail(errc::invalid_configuration);
    if (nearest < 0 || nearest >= static_cast<long long>(domain.cells[2]))
        return fail(errc::out_of_range);
    return static_cast<std::size_t>(nearest);
}

result<problem> make_problem(const mask_problem& value, std::optional<double> incident_power) noexcept {
    return detail::guard([&]() -> result<problem> {
        auto index = plane_index(value.grid, value.monitor_z);
        if (!index) return fail(index.error());

        problem out;
        out.grid = value.grid;
        out.materials = value.materials;
        out.geometry = value.mask;
        out.boundaries = value.boundaries;
        out.source = value.source;
        out.run = value.run;

        monitor detector;
        detector.name = "mask";
        detector.normal = axis::z;
        detector.index = *index;
        detector.time_window = value.monitor_time_window;
        auto options = value.spectrum;
        options.incident_power = incident_power;
        detector.spectrum = options;
        out.monitors.push_back(std::move(detector));
        return out;
    });
}

} // namespace

result<mask_spectrum> compute_mask_spectrum(const mask_problem& value) noexcept {
    return detail::guard([&]() -> result<mask_spectrum> {
        std::optional<spectrum_result> reference;
        double incident_power = 0;

        if (value.use_reference) {
            // First run the unpatterned stack to obtain the incident power and
            // the reference complex spectrum at the same monitor plane.
            auto reference_problem = make_problem(value, std::nullopt);
            if (!reference_problem) return fail(reference_problem.error());
            reference_problem->geometry.solids.clear();
            auto reference_run = solve(*reference_problem);
            if (!reference_run) return fail(reference_run.error());
            if (reference_run->monitors.empty() || !reference_run->monitors.front().spectrum)
                return fail(errc::solve_failed);
            reference = std::move(reference_run->monitors.front().spectrum);
            const auto power = zero_order_power(*reference);
            if (!power || !(*power > 0)) return fail(errc::numerical_failure);
            incident_power = *power;
        } else {
            incident_power = value.spectrum.incident_power.value_or(0);
        }

        auto full_problem = make_problem(value, incident_power > 0 ? std::optional<double>{incident_power}
                                                                  : std::nullopt);
        if (!full_problem) return fail(full_problem.error());
        auto full_run = solve(*full_problem);
        if (!full_run) return fail(full_run.error());
        if (full_run->monitors.empty() || !full_run->monitors.front().spectrum)
            return fail(errc::solve_failed);

        mask_spectrum out;
        out.spectrum = std::move(*full_run->monitors.front().spectrum);
        if (reference) {
            auto normalized = normalize_spectrum(out.spectrum, *reference);
            if (!normalized) return fail(normalized.error());
            out.normalized = std::move(*normalized);
        }
        out.reference = std::move(reference);
        out.incident_power = incident_power;
        out.total_efficiency = out.spectrum.total_efficiency.value_or(0);
        out.dt = full_run->dt;
        out.steps = full_run->steps;
        out.source_wavevector = full_run->source_wavevector;
        return out;
    });
}

} // namespace litho::fdtd
