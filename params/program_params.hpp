#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "particle_spec.hpp"

#include <array>

inline constexpr auto nThreads = 32;

inline constexpr auto x_collapsed = false;
inline constexpr auto y_collapsed = true;
inline constexpr auto z_collapsed = false;

inline constexpr auto Nx = 1501zu;
inline constexpr auto Ny = 2zu;
inline constexpr auto Nz = 1501zu;

inline constexpr std::array x_range = {-1.5e-05, 1.5e-05};
inline constexpr std::array y_range = {0.0, 0.01};
inline constexpr std::array z_range = {-1.5e-05, 1.5e-05};

inline constexpr auto dx = 2e-08;
inline constexpr auto dy = 0.01;
inline constexpr auto dz = 2e-08;

inline constexpr auto cfl   = 0.8479411200023808;
inline constexpr auto dt    = 4e-17;
inline constexpr auto t_end = 3e-13;
inline constexpr auto Nt    = 7500zu;

inline constexpr auto sim_name = "lsi_test";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = true;
inline constexpr auto push_enabled = true;
inline constexpr auto jdep_enabled = true;
inline constexpr auto coll_enabled = false;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 75zu;

inline constexpr auto PMLDepth    = 15zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {1zu, 1zu, 2zu, 2zu, 1zu, 1zu};

inline constexpr auto applied_fields_path = "";

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Reflecting, Periodic, Outflow };

inline constexpr auto particle_save_interval = 75zu;
inline constexpr auto interpolation_order = 1zu;

inline constexpr auto PBCSelect = ParticleBCType::Outflow;

inline constexpr std::array particle_spec = {
   ParticleGroupSpec{
      .name = "electrons",
      .filepath = "/data/lsi_test/electrons.bp",
      .mass = 9.1093837015e-31,
      .charge = -1.0,
      .atomic_number = 1
   },
   ParticleGroupSpec{
      .name = "ions",
      .filepath = "/data/lsi_test/ions.bp",
      .mass = 1.67262192369e-27,
      .charge = 1.0,
      .atomic_number = 1
   }
};

inline constexpr std::array<CollisionSpec, 0> collision_spec = {

};

/*---------------------------------------------------------------/
/-                      Metrics Parameters                      -/
/---------------------------------------------------------------*/
enum class MetricType { ParticleDump, ParticleDiag, ParticleEnergy, FieldDump, FieldEnergy };

inline constexpr auto metric_data_path = "/home/cepheid/TriForce/game_engine/data/lsi_test";
inline constexpr std::array<MetricType, 2> metric_spec = {
	MetricType::ParticleEnergy,
	MetricType::FieldEnergy
};

#endif //PROGRAM_PARAM_HPP
