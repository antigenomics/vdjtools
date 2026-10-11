#include "vdjtools/core.hpp"
#include "../../src/thread_joiner.hpp"
#include <stdexcept>

#include <cassert>
#include <cstdio>
#include <string>

// VDJTOOLS_VERSION is the single source of truth (parsed from pyproject.toml by CMake and
// injected into vdjtools_core as a PUBLIC compile definition). Asserting against the macro
// rather than a hand-copied literal verifies that version() actually returns the compiled-in
// version, and can never drift out of sync on a release bump — the failure this replaces.
int main() {
    using namespace vdjtools;
    assert(std::string(version()) == VDJTOOLS_VERSION);
    assert(std::string(version()).find('.') != std::string::npos);  // looks like X.Y.Z
    bool completed = false;
    try {
        std::vector<std::thread> pool;
        ThreadJoiner joiner{pool};
        pool.emplace_back([&] { completed = true; });
        throw std::runtime_error("later worker creation failed");
    } catch (const std::runtime_error&) {}
    assert(completed);  // no terminate, and shared state is alive until the worker finishes
    std::puts("ok");
    return 0;
}
