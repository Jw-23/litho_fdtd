#pragma once

#include <litho_fdtd/tracing.hpp>
#include <litho_fdtd/error.hpp>
#include <optional>

namespace litho::fdtd::detail {

// Logging is best effort at noexcept public boundaries, including allocation
// failures and synchronous exceptions from the tracing library.
class trace_span {
public:
    trace_span(trace_options options, const char* name, litho::Source source,
               const char* category = "fdtd") noexcept : options_(options), category_(category) {
        if (!options.session) return;
        try {
            if (options.session->category_enabled(category)) {
                scope_.emplace(*options.session, options.context, category, name, litho::Fields{}, source);
                options_.context = scope_->context();
            }
        } catch (...) {}
    }
    void finish() noexcept { scope_.reset(); }
    trace_options options() const noexcept { return options_; }
    template<class Fields>
    void event(const char* name, Fields&& fields, litho::Source source) noexcept {
        if (!options_.session) return;
        try {
            if (options_.session->category_enabled(category_))
                options_.session->emit(options_.context, category_, name, fields(), source);
        } catch (...) {}
    }
    void failure(errc error, litho::Source source) noexcept {
        if (scope_) scope_->fail();
        event("failed", [&] { return litho::Fields{
            litho::field("error_code", static_cast<int>(error)),
            litho::field("message", std::string(message(error)))}; }, source);
    }
private:
    trace_options options_;
    const char* category_;
    std::optional<litho::Scope> scope_;
};

} // namespace litho::fdtd::detail
