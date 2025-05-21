#ifndef TRIFORCE_TIMER_H
#define TRIFORCE_TIMER_H

#include <chrono>

//========== Timer Classes =========
//==================================
struct Timer
{
   using clock      = std::chrono::steady_clock;
   using duration_t = std::chrono::nanoseconds;
   using time_t     = std::chrono::time_point<clock>;

   time_t     start{};
   duration_t elapsed{};

   Timer()
   : start{clock::now()},
     elapsed{duration_t::zero()} {}

   void start_timer()
   {
      start = clock::now();
   }

   // For class timers that accumulate time spent class functions
   void stop_timer()
   {
      elapsed += clock::now() - start;
   }

   // For main timers that accumulates total runtime (t_start is never reset)
   void update_duration()
   {
      elapsed = clock::now() - start;
   }

   void reset()
   {
      start   = clock::now();
      elapsed = duration_t::zero();
   }
};

#endif //TRIFORCE_TIMER_H
