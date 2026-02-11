#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "particle_spec.hpp"

#include <array>

inline constexpr auto nThreads = 48;

inline constexpr auto x_collapsed = false;
inline constexpr auto y_collapsed = false;
inline constexpr auto z_collapsed = false;

inline constexpr auto Nx = 57zu;
inline constexpr auto Ny = 57zu;
inline constexpr auto Nz = 113zu;

inline constexpr std::array x_range = {-0.14, 0.14};
inline constexpr std::array y_range = {-0.14, 0.14};
inline constexpr std::array z_range = {-0.28, 0.28};

inline constexpr auto dx = 0.005;
inline constexpr auto dy = 0.005;
inline constexpr auto dz = 0.005;

inline constexpr auto cfl   = 0.5192557689819587;
inline constexpr auto dt    = 5e-12;
inline constexpr auto t_end = 8e-08;
inline constexpr auto Nt    = 16001zu;

inline constexpr auto sim_name = "seinfeld3D";
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

inline constexpr auto em_save_interval = 20zu;

inline constexpr auto PMLDepth    = 6zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {1zu, 1zu, 1zu, 1zu, 1zu, 1zu};

inline constexpr auto applied_fields_path = "/home/cepheid/TriForce/game_engine/data/seinfeld3D_applied_fields.bp";

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Reflecting, Periodic, Outflow };

inline constexpr auto particle_save_interval = 20zu;
inline constexpr auto interpolation_order = 1zu;

inline constexpr auto PBCSelect = ParticleBCType::Outflow;

inline constexpr std::array particle_spec = {
   ParticleGroupSpec{
      .name = "electrons",
      .filepath = "/data/electrons.bp",
      .mass = 9.1093837015e-31,
      .charge = -1.0,
      .atomic_number = 0
   },
   ParticleGroupSpec{
      .name = "ions",
      .filepath = "/data/ions.bp",
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

inline constexpr auto metric_data_path = "/home/cepheid/TriForce/game_engine/data";
inline constexpr std::array<MetricType, 2> metric_spec = {
	MetricType::ParticleEnergy,
	MetricType::FieldEnergy
};

#endif //PROGRAM_PARAM_HPP
