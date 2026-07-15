#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "particle_spec.hpp"
#include "sources_spec.hpp"

#include <array>

inline constexpr auto nThreads = 64;

inline constexpr auto x_collapsed = false;
inline constexpr auto y_collapsed = true;
inline constexpr auto z_collapsed = false;

inline constexpr auto Nx = 1551zu;
inline constexpr auto Ny = 2zu;
inline constexpr auto Nz = 1551zu;

inline constexpr std::array x_range = {-1.55e-05, 1.55e-05};
inline constexpr std::array y_range = {0.0, 0.01};
inline constexpr std::array z_range = {-1.55e-05, 1.55e-05};

inline constexpr auto dx = 2e-08;
inline constexpr auto dy = 0.01;
inline constexpr auto dz = 2e-08;

inline constexpr auto dt    = 4e-17;
inline constexpr auto t_end = 3e-13;
inline constexpr auto Nt    = 7500zu;

inline constexpr auto sim_name = "lsi_smith";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = true;
inline constexpr auto push_enabled = true;
inline constexpr auto jdep_enabled = true;
inline constexpr auto coll_enabled = false;
inline constexpr auto applied_fields_only = false;
inline constexpr auto velocity_backstep_enabled = true;
inline constexpr auto ionization_test_enabled = false;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 75zu;
inline constexpr auto em_subcycles = 1zu;
inline constexpr auto dt_em = 4e-17;

inline constexpr auto PMLDepth    = 15zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {1zu, 1zu, 2zu, 2zu, 1zu, 1zu};

inline constexpr auto laser_enabled = true;
inline constexpr auto applied_fields_path = "";

constexpr auto laser_spec = LaserSpec{.lambda=8e-07, .E0=-27500000000000.0, .w0=2.5479e-06, .xspot=1.5e-05, .scale=1.28855495};

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Reflecting, Periodic, Outflow };
enum class ParticlePushType { Ballistic, Boris, HigueraCary };

inline constexpr auto particle_save_interval = 75zu;
inline constexpr auto sort_frequency = 100zu;
inline constexpr auto interpolation_order = 1zu;
inline constexpr auto ParticlePushSelect = ParticlePushType::Boris;
inline constexpr auto PBCSelect = ParticleBCType::Outflow;
inline constexpr auto PBCDepth = 3zu;

inline constexpr std::array<ParticleGroupSpec, 2> particle_spec = {
   ParticleGroupSpec{
      .name = "hydrogen",
      .filepath = "/data/lsi_smith/hydrogen.bp",
      .mass = 1.67382338147136e-27,
      .charge = 1.0,
      .atomic_number = 1,
      .tracer = false
   },
   ParticleGroupSpec{
      .name = "electrons",
      .filepath = "/data/lsi_smith/electrons.bp",
      .mass = 9.1093837139e-31,
      .charge = -1.0,
      .atomic_number = 0,
      .tracer = false
   }
};

inline constexpr std::array<CollisionSpec, 0> collision_spec = {

};

/*---------------------------------------------------------------/
/-                      Metrics Parameters                      -/
/---------------------------------------------------------------*/
enum class MetricType { ParticleDump, ParticleDiag, ParticleEnergy, FieldDump, FieldEnergy };

inline constexpr auto metric_data_path = "/home/cepheid/TriForce/game_engine/data/lsi_smith";
inline constexpr std::array<MetricType, 3> metric_spec = {
	MetricType::ParticleEnergy,
	MetricType::FieldEnergy,
	MetricType::ParticleDiag
};

#endif //PROGRAM_PARAM_HPP
