#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "particle_spec.hpp"

#include <array>

inline constexpr auto nThreads = 4;

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

inline constexpr auto cfl   = 0.2596278844909794;
inline constexpr auto dt    = 5e-16;
inline constexpr auto t_end = 5e-12;
inline constexpr auto Nt    = 10000zu;

inline constexpr auto sim_name = "carbon_thermal_eq";
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

inline constexpr auto particle_save_interval = 100zu;
inline constexpr auto interpolation_order = 1zu;

inline constexpr auto PBCSelect = ParticleBCType::Periodic;

inline constexpr std::array particle_spec = {
   ParticleGroupSpec{
      .name = "carbon1",
      .filepath = "/data/carbon1.bp",
      .mass = 1.9945e-26,
      .charge = 6.0,
      .atomic_number = 6
   },
   ParticleGroupSpec{
      .name = "carbon2",
      .filepath = "/data/carbon2.bp",
      .mass = 1.9945e-26,
      .charge = 6.0,
      .atomic_number = 6
   }
};

inline constexpr std::array collision_spec = {
   CollisionSpec{
      .group1 = "carbon1",
      .group2 = "carbon2",
      .channels = {"coulomb"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = false,
      .coulomb = {.coulomb_log = 10.0, .rate_multiplier = 1.0},						
   },
   CollisionSpec{
      .group1 = "carbon1",
      .group2 = "carbon1",
      .channels = {"coulomb"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = true,
      .coulomb = {.coulomb_log = 10.0, .rate_multiplier = 1.0},						
   },
   CollisionSpec{
      .group1 = "carbon2",
      .group2 = "carbon2",
      .channels = {"coulomb"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = true,
      .coulomb = {.coulomb_log = 10.0, .rate_multiplier = 1.0},						
   }
};

#endif //PROGRAM_PARAM_HPP
