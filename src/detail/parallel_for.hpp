#pragma once

#include <litho_fdtd/error.hpp>

#include <parallel.hpp>

#include <algorithm>
#include <cstddef>
#include <system_error>
#include <thread>
#include <utility>

namespace litho::fdtd::detail {

[[nodiscard]] inline std::size_t worker_count(std::size_t workers) noexcept {
    if (workers != 0) return workers;
    return std::max<std::size_t>(1, std::thread::hardware_concurrency());
}

// Dispatch a cache-friendly chunk body through the project's parallel
// scheduler. Tiny chunk counts stay inline to avoid waking a team for work
// that is shorter than the synchronization cost. The body writes disjoint
// slices and does not allocate.
template<class Body>
[[nodiscard]] result<void> parallel_for(std::size_t chunks, std::size_t workers, Body&& body) {
    if (chunks == 0) return {};
    // An explicit worker count is honored even for small grids; automatic mode
    // stays inline when the team would cost more than the work.
    const std::size_t actual = workers == 1 ? 1 : (workers == 0 && chunks < 4096) ? 1 : worker_count(workers);
    auto status = parallel::for_each_chunk(chunks, actual, std::forward<Body>(body));
    if (!status) {
        if (status.error() == std::errc::not_enough_memory) return fail(errc::allocation_failure);
        return fail(errc::solve_failed);
    }
    return {};
}

} // namespace litho::fdtd::detail
