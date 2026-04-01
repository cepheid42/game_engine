#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "particle_spec.hpp"

#include <array>

inline constexpr auto nThreads = 1;

inline constexpr auto x_collapsed = false;
inline constexpr auto y_collapsed = false;
inline constexpr auto z_collapsed = false;

inline constexpr auto Nx = 512zu;
inline constexpr auto Ny = 512zu;
inline constexpr auto Nz = 128zu;

inline constexpr std::array x_range = {-200000000.0, 200000000.0};
inline constexpr std::array y_range = {-200000000.0, 200000000.0};
inline constexpr std::array z_range = {-20000000000000.0, 20000000000000.0};

inline constexpr auto dx = 782778.8649706458;
inline constexpr auto dy = 782778.8649706458;
inline constexpr auto dz = 314960629921.2598;

inline constexpr auto dt    = 1e-07;
inline constexpr auto t_end = 0.3141592653589793;
inline constexpr auto Nt    = 3141593zu;

inline constexpr auto sim_name = "magnetic_mirror_HigueraCary";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = false;
inline constexpr auto push_enabled = true;
inline constexpr auto jdep_enabled = false;
inline constexpr auto coll_enabled = false;
inline constexpr auto applied_fields_only = true;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 500000zu;
inline constexpr auto em_subcycles = 1zu;
inline constexpr auto dt_em = 0.0018278417169303963;

inline constexpr auto PMLDepth    = 10zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {2zu, 2zu, 2zu, 2zu, 2zu, 2zu};

inline constexpr auto applied_fields_path = "/home/cepheid/TriForce/game_engine/data/magnetic_mirror_HigueraCary/magnetic_mirror_HigueraCary_applied_fields.bp";

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Reflecting, Periodic, Outflow };
enum class ParticlePushType { Ballistic, Boris, HigueraCary };

inline constexpr auto particle_save_interval = 500000zu;
inline constexpr auto sort_frequency = 100zu;
inline constexpr auto interpolation_order = 1zu;
inline constexpr auto ParticlePushSelect = ParticlePushType::HigueraCary;
inline constexpr auto PBCSelect = ParticleBCType::Outflow;
inline constexpr auto PBCDepth = 0zu;

inline constexpr std::array<ParticleGroupSpec, 1> particle_spec = {
   ParticleGroupSpec{
      .name = "sp",
      .filepath = "/data/magnetic_mirror_HigueraCary/sp.bp",
      .mass = 1.0,
      .charge = 6.241509074460763e+18,
      .atomic_number = 0,
      .tracer = true
   }
};

inline constexpr std::array<CollisionSpec, 0> collision_spec = {

};

/*---------------------------------------------------------------/
/-                      Metrics Parameters                      -/
/---------------------------------------------------------------*/
enum class MetricType { ParticleDump, ParticleDiag, ParticleEnergy, FieldDump, FieldEnergy };

inline constexpr auto metric_data_path = "/home/cepheid/TriForce/game_engine/data/magnetic_mirror_HigueraCary";
inline constexpr std::array<MetricType, 1> metric_spec = {
	MetricType::ParticleDump
};

#endif //PROGRAM_PARAM_HPP
