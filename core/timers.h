//
// Created by cepheid on 11/9/23.
//

#ifndef TRIFORCE_TIMER_H
#define TRIFORCE_TIMER_H

#include <chrono>
#include <map>
#include <string>
#include <iostream>

//========== Timer Classes =========
//==================================
struct Timer {
  using clock = std::chrono::steady_clock;
  using duration_t = std::chrono::nanoseconds;
  using time_t = std::chrono::time_point<clock>;

  duration_t elapsed{};
  time_t start{};

  Timer()
  : elapsed{duration_t::zero()},
    start{clock::now()}
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
  }
};

struct EventTimer {
  Timer& timer;

  explicit EventTimer(Timer& timer_)
  : timer{timer_}
  {
    timer.start_timer();
  }

  ~EventTimer() {
    timer.stop_timer();
  }
}

#endif //TRIFORCE_TIMER_H
