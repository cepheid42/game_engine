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
struct CoulombSpec {
   double coulomb_log{10.0};
   double rate_multiplier{1.0};
};

struct IonizationSpec {
   std::string_view product1{};
   std::string_view product2{}; // no product 3 for now
   std::string_view cross_section_file{};
   double ionization_energy{0.0};
   double rate_multiplier{1.0};
   double production_multiplier{1.0};
   double rejection_multiplier{1.0};
   double constant_cross_section{0.0};
};

struct FusionSpec {
   std::string_view product1{};
   std::string_view product2{};
   std::string_view cross_section_file{};
   double energy_gain{0.0};
   double rate_multiplier{1.0};
   double production_multiplier{1.0};
   double constant_cross_section{0.0};
};

struct RadiationSpec {
   std::string_view product1{};
   std::string_view cross_section_file{};
   double production_multiplier{1.0};
   double rate_multiplier{1.0};
   bool reduce_electron_energy{false};
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
   RadiationSpec radiation{};
};

#endif //GAME_ENGINE_PARTICLE_SPEC_HPP