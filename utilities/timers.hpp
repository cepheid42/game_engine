#ifndef TRIFORCE_TIMER_H
#define TRIFORCE_TIMER_H

#include <chrono>
#include <print>
#include <unordered_map>
#include <string>

namespace tf::utilities {
//========== Timer Classes =========
//==================================
struct Timer
{
   using clock      = std::chrono::steady_clock;
   using duration_t = std::chrono::nanoseconds;
   using time_t     = std::chrono::time_point<clock>;

   time_t     start{clock::now()};
   duration_t elapsed{duration_t::zero()};

   void start_timer() { start = clock::now(); }

   // For timers that accumulate time spent in functions
   void stop_timer() { elapsed += clock::now() - start; }

   // For timers that accumulate total runtime (t_start is never reset)
   void update_duration() { elapsed = clock::now() - start;}

   void reset() {
      start   = clock::now();
      elapsed = duration_t::zero();
   }

   [[nodiscard]] auto to_seconds() const {
      return std::chrono::duration<double>{elapsed};
   }
}; // end struct Timer

inline auto create_timers() -> std::unordered_map<std::string, Timer> {
   std::unordered_map<std::string, Timer> timers{};
   timers["Main"] = Timer{};
   timers["EM"] = Timer{};
   timers["Push"] = Timer{};
   timers["Jdep"] = Timer{};
   timers["Collisions"] = Timer{};
   timers["IO"] = Timer{};
   return timers;
} // end create_timers()

inline void print_final_timers(std::unordered_map<std::string, Timer>& timers) {
   std::println("        EM: {}", std::chrono::hh_mm_ss(timers["EM"].elapsed));
   std::println("      Push: {}", std::chrono::hh_mm_ss(timers["Push"].elapsed));
   std::println("      Jdep: {}", std::chrono::hh_mm_ss(timers["Jdep"].elapsed));
   std::println("        IO: {}", std::chrono::hh_mm_ss(timers["IO"].elapsed));
   std::println("Collisions: {}", std::chrono::hh_mm_ss(timers["Collisions"].elapsed));
   std::println("     Total: {}", std::chrono::hh_mm_ss(timers["Main"].elapsed));
}

} // end namespace tf::utilities
#endif //TRIFORCE_TIMER_H
