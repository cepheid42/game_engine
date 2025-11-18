#ifndef GAME_ENGINE_PARTICLE_PARAMS_HPP
#define GAME_ENGINE_PARTICLE_PARAMS_HPP

#include <string>

struct CollisionParams {
   consteval CollisionParams() = delete;
   constexpr CollisionParams(const std::string_view g1_name_, const std::string_view g2_name_, const std::size_t interval_, const double clog_, const double rate_)
   : group1_name(g1_name_), group2_name(g2_name_), step_interval(interval_), coulombLog(clog_), rate_mult(rate_)
   {}

   std::string_view group1_name;
   std::string_view group2_name;
   std::size_t step_interval;
   double coulombLog;
   double rate_mult;
};


#endif //GAME_ENGINE_PARTICLE_PARAMS_HPP