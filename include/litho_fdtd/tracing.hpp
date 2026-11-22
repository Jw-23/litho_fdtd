#pragma once

#include <litho_trace/trace.hpp>
#include <cstddef>

namespace litho::fdtd {

// Non-owning, explicit trace context. The session must outlive the solve.
// A null session disables instrumentation. Trace failures never fail a solve;
// inspect Session::stats() after closing the session for dropped events/I/O errors.
struct trace_options {
    litho::Session* session = nullptr;
    litho::Context context{};
    std::size_t step_interval = 100; // 0 disables step spans
};

} // namespace litho::fdtd
