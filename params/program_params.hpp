#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "particle_spec.hpp"

#include <array>

inline constexpr auto nThreads = 1;

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

inline constexpr auto cfl   = 25.962788449097935;
inline constexpr auto dt    = 5e-14;
inline constexpr auto t_end = 1e-12;
inline constexpr auto Nt    = 21zu;

inline constexpr auto sim_name = "DT_fusion_5keV";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = false;
inline constexpr auto push_enabled = false;
inline constexpr auto jdep_enabled = false;
inline constexpr auto coll_enabled = true;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 1zu;

inline constexpr auto PMLDepth    = 10zu;
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

inline constexpr std::array particle_spec = {
   ParticleGroupSpec{
      .name = "tritium",
      .filepath = "/data/tritium.bp",
      .mass = 5.0082688518189296e-27,
      .charge = 1.0,
      .atomic_number = 1
   },
   ParticleGroupSpec{
      .name = "deuterium",
      .filepath = "/data/deuterium.bp",
      .mass = 3.3234028878932395e-27,
      .charge = 1.0,
      .atomic_number = 1
   },
   ParticleGroupSpec{
      .name = "neutrons",
      .filepath = "",
      .mass = 1.67492749804e-27,
      .charge = 0.0,
      .atomic_number = 0
   },
   ParticleGroupSpec{
      .name = "helium4",
      .filepath = "",
      .mass = 6.64647864959036e-27,
      .charge = 2.0,
      .atomic_number = 2
   }
};

inline constexpr std::array collision_spec = {
   CollisionSpec{
      .group1 = "deuterium",
      .group2 = "tritium",
      .channels = {"fusion"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = false,
      .fusion = {
         .product1 = "neutrons",
         .product2 = "helium4",
         .energy_gain = 17589000.0,
         .rate_multiplier = 1.0,
         .production_multiplier = 10000000000.0,
         .constant_cross_section = 0.0,
         .cross_section_file = "/data/DT_nHe4_BH_eV_m2.txt",
      },
   }
};

#endif //PROGRAM_PARAM_HPP
