#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "particle_spec.hpp"

#include <array>

inline constexpr auto nThreads = 32;

inline constexpr auto x_collapsed = false;
inline constexpr auto y_collapsed = false;
inline constexpr auto z_collapsed = true;

inline constexpr auto Nx = 213zu;
inline constexpr auto Ny = 213zu;
inline constexpr auto Nz = 2zu;

inline constexpr std::array x_range = {-0.03, 1.03};
inline constexpr std::array y_range = {-0.03, 1.03};
inline constexpr std::array z_range = {0.0, 1.0};

inline constexpr auto dx = 0.005;
inline constexpr auto dy = 0.005;
inline constexpr auto dz = 1.0;

inline constexpr auto cfl   = 211.98660490424294;
inline constexpr auto dt    = 2.5e-09;
inline constexpr auto t_end = 5e-06;
inline constexpr auto Nt    = 2001zu;

inline constexpr auto sim_name = "efield_only";
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

inline constexpr auto em_save_interval = 50zu;
inline constexpr auto em_subcycles = 215zu;
inline constexpr auto dt_em = 1.1675265996726487e-11;

inline constexpr auto PMLDepth    = 6zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {1zu, 1zu, 1zu, 1zu, 2zu, 2zu};

inline constexpr auto applied_fields_path = "";

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Reflecting, Periodic, Outflow };

inline constexpr auto particle_save_interval = 50zu;
inline constexpr auto interpolation_order = 1zu;

inline constexpr auto PBCSelect = ParticleBCType::Outflow;

inline constexpr std::array<ParticleGroupSpec, 2> particle_spec = {
   ParticleGroupSpec{
      .name = "electrons",
      .filepath = "/data/efield_only/electrons.bp",
      .mass = 9.1093837015e-31,
      .charge = -1.0,
      .atomic_number = 0
   },
   ParticleGroupSpec{
      .name = "ions",
      .filepath = "/data/efield_only/ions.bp",
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

inline constexpr auto metric_data_path = "/home/cepheid/TriForce/game_engine/data/efield_only";
inline constexpr std::array<MetricType, 2> metric_spec = {
	MetricType::ParticleDump,
	MetricType::FieldEnergy
};

#endif //PROGRAM_PARAM_HPP
