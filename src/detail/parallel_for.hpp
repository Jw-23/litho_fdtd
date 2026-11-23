#pragma once

#include "tracing.hpp"

#include <litho_fdtd/error.hpp>

#include <parallel.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <system_error>
#include <thread>
#include <utility>

namespace litho::fdtd::detail {

// Match the development M4's four performance cores with a four-worker budget.
// This is a thread-count limit, not CPU affinity: the OS still chooses cores.
inline constexpr std::size_t max_workers = 4;

[[nodiscard]] inline std::size_t worker_count(std::size_t workers) noexcept {
    if (workers == 0) workers = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    return std::min(max_workers, std::max<std::size_t>(1, workers));
}

[[nodiscard]] inline std::size_t effective_workers(std::size_t chunks, std::size_t workers) noexcept {
    return workers == 1 || (workers == 0 && chunks < 4096) ? 1 : worker_count(workers);
}

// Dispatch a cache-friendly chunk body through the project's parallel
// scheduler. Tiny chunk counts stay inline to avoid waking a team for work
// that is shorter than the synchronization cost. The body writes disjoint
// slices and does not allocate.
//
// When a trace session is attached, the kernel span carries a summary event:
// wall_ns, main_work_ns, wait_ns, total work_ns, chunk count and worker
// utilization. The wait field is the calling thread's time outside its own
// chunks; utilization distinguishes barrier/idle time from memory contention.
// Chunk spans are emitted only for the first trace.detail_steps sampled steps.
template<class Body>
[[nodiscard]] result<void> parallel_for(std::size_t work_items, std::size_t workers, std::size_t chunk_rows,
                                        const char* phase, trace_options trace, std::size_t step, Body&& body) {
    if (work_items == 0) return {};
    const std::size_t actual = effective_workers(work_items, workers);
    auto map_error = [](const parallel::Result<void>& status) -> result<void> {
        if (!status) {
            if (status.error() == std::errc::not_enough_memory) return fail(errc::allocation_failure);
            return fail(errc::solve_failed);
        }
        return {};
    };
    if (!trace.session)
        return map_error(parallel::for_each_chunk(work_items, actual, std::forward<Body>(body), chunk_rows));

    detail::trace_span phase_span(trace, phase, LT_SOURCE, "fdtd.kernel");
    const bool detailed = step < trace.detail_steps;
    const auto main_thread = std::this_thread::get_id();
    const auto start = std::chrono::steady_clock::now();
    std::uint64_t main_work_ns = 0;
    std::atomic<std::uint64_t> total_work_ns{0};
    std::atomic<std::size_t> chunk_count{0};

    auto status = parallel::for_each_chunk(work_items, actual, [&](std::size_t first, std::size_t last) {
        if (detailed) {
            detail::trace_span chunk_span(phase_span.options(), "chunk", LT_SOURCE, "fdtd.chunk");
            const auto chunk_start = std::chrono::steady_clock::now();
            body(first, last);
            const auto nanoseconds = static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() -
                                                                     chunk_start)
                    .count());
            total_work_ns.fetch_add(nanoseconds, std::memory_order_relaxed);
            chunk_count.fetch_add(1, std::memory_order_relaxed);
            if (std::this_thread::get_id() == main_thread) main_work_ns += nanoseconds;
        } else if (std::this_thread::get_id() == main_thread) {
            const auto chunk_start = std::chrono::steady_clock::now();
            body(first, last);
            main_work_ns += static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() -
                                                                     chunk_start)
                    .count());
        } else {
            body(first, last);
        }
    }, chunk_rows);

    const auto wall_ns = static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count());
    const std::uint64_t wait_ns = wall_ns > main_work_ns ? wall_ns - main_work_ns : 0;
    const double utilization = actual != 0 && wall_ns != 0
                                   ? static_cast<double>(total_work_ns.load(std::memory_order_relaxed)) /
                                         (static_cast<double>(wall_ns) * static_cast<double>(actual))
                                   : 0.0;
    const double main_fraction =
        wall_ns != 0 ? static_cast<double>(main_work_ns) / static_cast<double>(wall_ns) : 0.0;
    phase_span.event(
        "summary",
        [&] {
            return litho::Fields{litho::field("step", static_cast<std::uint64_t>(step)),
                                 litho::field("workers", static_cast<std::uint64_t>(actual)),
                                 litho::field("items", static_cast<std::uint64_t>(work_items)),
                                 litho::field("wall_ns", wall_ns),
                                 litho::field("main_work_ns", main_work_ns),
                                 litho::field("wait_ns", wait_ns),
                                 litho::field("chunks", static_cast<std::uint64_t>(chunk_count.load(
                                                           std::memory_order_relaxed))),
                                 litho::field("work_ns", total_work_ns.load(std::memory_order_relaxed)),
                                 litho::field("utilization", utilization),
                                 litho::field("main_fraction", main_fraction)};
        },
        LT_SOURCE);
    return map_error(status);
}

} // namespace litho::fdtd::detail
