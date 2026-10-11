#pragma once
#include <thread>
#include <vector>

namespace vdjtools {
// Join already-started workers if a later thread creation throws, before shared state unwinds.
struct ThreadJoiner {
    std::vector<std::thread>& threads;
    ~ThreadJoiner() {
        for (auto& thread : threads) if (thread.joinable()) thread.join();
    }
};
}  // namespace vdjtools
