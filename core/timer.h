//
// Created by cepheid on 7/15/24.
//

#ifndef TIMER_H
#define TIMER_H

#include <chrono>

struct EventTimer {
    using clock = std::chrono::steady_clock;
    using duration_t = std::chrono::nanoseconds;
    using time_t = std::chrono::time_point<clock>;

    time_t start{};
    duration_t elapsed{};

    EventTimer()
    : start{clock::now()},
      elapsed{duration_t::zero()}
    {}

    void start_timer() {
        start = clock::now();
    }

    // For class timers that accumulate time spent class functions
    void stop_timer() {
        elapsed += clock::now() - start;
    }

    // For main timers that accumulates total runtime (t_start is never reset)
    void update_duration() {
        elapsed = clock::now() - start;
    }

    void reset() {
        start     = clock::now();
        elapsed   = duration_t::zero();
    };
};

#endif //TIMER_H
