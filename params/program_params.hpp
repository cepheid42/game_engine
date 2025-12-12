#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>
#include <string>
#include <tuple>

inline constexpr auto nThreads = 8;

inline constexpr auto x_collapsed = true;
inline constexpr auto y_collapsed = true;
inline constexpr auto z_collapsed = true;

inline constexpr auto Nx = 2zu;
inline constexpr auto Ny = 2zu;
inline constexpr auto Nz = 2zu;

inline constexpr std::array x_range = {0.0, 1e-06};
inline constexpr std::array y_range = {0.0, 1e-06};
inline constexpr std::array z_range = {0.0, 1e-06};

inline constexpr auto dx = 1e-06;
inline constexpr auto dy = 1e-06;
inline constexpr auto dz = 1e-06;

inline constexpr auto cfl   = 0.0025962788449097936;
inline constexpr auto dt    = 5e-18;
inline constexpr auto t_end = 3.18e-15;
inline constexpr auto Nt    = 637zu;

inline constexpr auto sim_name = "ionization";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = false;
inline constexpr auto push_enabled = true;
inline constexpr auto jdep_enabled = false;
inline constexpr auto coll_enabled = true;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 1zu;

inline constexpr auto PMLDepth    = 0zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {2zu, 2zu, 2zu, 2zu, 2zu, 2zu};

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Reflecting, Periodic, Outflow };

inline constexpr auto particle_save_interval = 1zu;
inline constexpr auto interpolation_order = 1zu;

inline constexpr auto PBCSelect = ParticleBCType::Periodic;

struct ParticleGroupSpec {
   std::string_view name;
   std::string_view filepath;
   double mass;
   double charge;
   std::size_t atomic_number;
};

inline constexpr std::array particle_spec = {
   ParticleGroupSpec{"electrons", "/data/electrons.bp", 9.1093837015e-31, -1.602176634e-19, 0},
   ParticleGroupSpec{"Al", "/data/Al.bp", 4.481567702427985e-26, 0.0, 13},
   ParticleGroupSpec{"Al+", "", 4.4814766085909697e-26, 1.602176634e-19, 13}
};

struct IonizationSpec {
   std::string product1{};
   std::string product2{}; // no product 3 for now
   double ionization_energy{1.0};
   double rate_multiplier{1.0};
   double production_multiplier{1.0};
   double rejection_multiplier{1.0};
   double constant_cross_section{0.0};
   std::string cross_section_file{};
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
};

inline constexpr std::array collision_spec = {
   CollisionSpec{
      .group1 = "electrons",
      .group2 = "Al",
      .channels = {"ionization"}, // no coulombs in this example
      .step_interval = 1zu,
      .probability_search_area = 1.0,
      .self_scatter = false,
      // .coulomb = {.coulomb_log = 10.0, .rate_multiplier = 1.0},
      .ionization = {
         .product1 = "electrons",
         .product2 = "Al+",
         .ionization_energy = 5.9858,
         .rate_multiplier = 1.0,
         .production_multiplier = 1.0,
         .rejection_multiplier = 1.0,
         .constant_cross_section = 0.0,
         .cross_section_file = "",
      }
   }
};

#endif //PROGRAM_PARAM_HPP
