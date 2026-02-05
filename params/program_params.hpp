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

inline constexpr auto cfl   = 0.5192557689819588;
inline constexpr auto dt    = 1e-15;
inline constexpr auto t_end = 4.2e-14;
inline constexpr auto Nt    = 42zu;

inline constexpr auto sim_name = "2wiwe";
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
      .name = "electrons",
      .filepath = "/data/electrons.bp",
      .mass = 9.1093837015e-31,
      .charge = -1.0,
      .atomic_number = 0
   },
   ParticleGroupSpec{
      .name = "copper",
      .filepath = "/data/copper.bp",
      .mass = 1.0552725768242999e-25,
      .charge = 0.0,
      .atomic_number = 29
   },
   ParticleGroupSpec{
      .name = "photons",
      .filepath = "",
      .mass = 0.0,
      .charge = 0.0,
      .atomic_number = 0
   }
};

inline constexpr std::array collision_spec = {
   CollisionSpec{
      .group1 = "electrons",
      .group2 = "copper",
      .channels = {"radiation"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = false,
      .radiation = {
         .product1 = "photons",
         .cross_section_file = "/tests/cross_section_data/SB_G4_Z29_kdsdk_MeV_barns.csv",
         .production_multiplier = 100000.0,
         .reduce_electron_energy = false,
      },
   }
};

#endif //PROGRAM_PARAM_HPP
