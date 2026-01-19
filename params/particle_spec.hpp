#ifndef GAME_ENGINE_PARTICLE_SPEC_HPP
#define GAME_ENGINE_PARTICLE_SPEC_HPP

#include <array>
#include <string>

struct ParticleGroupSpec {
   std::string_view name;
   std::string_view filepath;
   double mass;
   double charge;
   std::size_t atomic_number;
};

struct FusionSpec {
   std::string_view product1{};
   std::string_view product2{};
   double energy_gain{0.0};
   double rate_multiplier{1.0};
   double production_multiplier{1.0};
   double constant_cross_section{0.0};
   std::string_view cross_section_file{};
};

struct IonizationSpec {
   std::string_view product1{};
   std::string_view product2{}; // no product 3 for now
   double ionization_energy{0.0};
   double rate_multiplier{1.0};
   double production_multiplier{1.0};
   double rejection_multiplier{1.0};
   double constant_cross_section{0.0};
   std::string_view cross_section_file{};
};

struct CoulombSpec {
   double coulomb_log{10.0};
   double rate_multiplier{1.0};
};

struct CollisionSpec {
   std::string_view group1;
   std::string_view group2;
   std::array<std::string_view, 2> channels;
   std::size_t step_interval;
   double probability_search_area;
   bool self_scatter;
   CoulombSpec coulomb{};
   IonizationSpec ionization{};
   FusionSpec fusion{};
};

#endif //GAME_ENGINE_PARTICLE_SPEC_HPP