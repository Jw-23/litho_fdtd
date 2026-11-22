#include <litho_fdtd/mask.hpp>
#include "detail/tracing.hpp"

#include <cmath>
#include <utility>

namespace litho::fdtd {

namespace {

result<problem> make_problem(const mask_problem& value, std::optional<double> incident_power) noexcept {
    return detail::guard([&]() -> result<problem> {
        auto plane = monitor_at(value.grid, axis::z, value.monitor_z);
        if (!plane) return fail(plane.error());

        problem out;
        out.grid = value.grid;
        out.materials = value.materials;
        out.geometry = value.mask;
        out.boundaries = value.boundaries;
        out.source = value.source;
        out.run = value.run;

        monitor detector = std::move(*plane);
        detector.name = "mask";
        detector.normal = axis::z;
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
    detail::trace_span span(value.run.trace, "mask_spectrum", LT_SOURCE);
    auto outcome = detail::guard([&]() -> result<mask_spectrum> {
        std::optional<spectrum_result> reference;
        double incident_power = 0;

        if (value.use_reference) {
            // First run the unpatterned stack to obtain the incident power and
            // the reference complex spectrum at the same monitor plane.
            auto reference_problem = make_problem(value, std::nullopt);
            if (!reference_problem) return fail(reference_problem.error());
            reference_problem->geometry.solids.clear();
            detail::trace_span reference_span(span.options(), "reference", LT_SOURCE);
            reference_problem->run.trace = reference_span.options();
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
        detail::trace_span full_span(span.options(), "patterned", LT_SOURCE);
        full_problem->run.trace = full_span.options();
        auto full_run = solve(*full_problem);
        if (!full_run) return fail(full_run.error());
        full_span.finish();
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
    if (!outcome) span.failure(outcome.error(), LT_SOURCE);
    return outcome;
}

} // namespace litho::fdtd
