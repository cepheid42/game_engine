#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "particle_spec.hpp"
#include "sources_spec.hpp"

#include <array>

inline constexpr auto nThreads = 48;

inline constexpr auto x_collapsed = false;
inline constexpr auto y_collapsed = true;
inline constexpr auto z_collapsed = false;

inline constexpr auto Nx = 551zu;
inline constexpr auto Ny = 2zu;
inline constexpr auto Nz = 151zu;

inline constexpr std::array x_range = {-1.55e-05, -4.5e-06};
inline constexpr std::array y_range = {0.0, 0.01};
inline constexpr std::array z_range = {-1.55e-05, 1.55e-05};

inline constexpr auto dx = 2e-08;
inline constexpr auto dy = 0.01;
inline constexpr auto dz = 2.0666666666666666e-07;

inline constexpr auto dt    = 2e-17;
inline constexpr auto t_end = 4e-17;
inline constexpr auto Nt    = 3zu;

inline constexpr auto sim_name = "rlsi_fusion";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = true;
inline constexpr auto push_enabled = true;
inline constexpr auto jdep_enabled = true;
inline constexpr auto coll_enabled = true;
inline constexpr auto applied_fields_only = false;
inline constexpr auto velocity_backstep_enabled = true;
inline constexpr auto ionization_test_enabled = false;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 1zu;
inline constexpr auto em_subcycles = 1zu;
inline constexpr auto dt_em = 2e-17;

inline constexpr auto PMLDepth    = 15zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {1zu, 1zu, 2zu, 2zu, 1zu, 1zu};

inline constexpr auto laser_enabled = true;
inline constexpr auto applied_fields_path = "";

constexpr auto laser_spec = LaserSpec{.lambda=8e-07, .E0=-27500000000000.0, .w0=2.5479e-06, .xspot=5e-06, .scale=0.60454};

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Reflecting, Periodic, Outflow };
enum class ParticlePushType { Ballistic, Boris, HigueraCary };

inline constexpr auto particle_save_interval = 1zu;
inline constexpr auto sort_frequency = 100zu;
inline constexpr auto interpolation_order = 1zu;
inline constexpr auto ParticlePushSelect = ParticlePushType::Boris;
inline constexpr auto PBCSelect = ParticleBCType::Outflow;
inline constexpr auto PBCDepth = 3zu;

inline constexpr std::array<ParticleGroupSpec, 7> particle_spec = {
   ParticleGroupSpec{
      .name = "deuterium",
      .filepath = "/data/rlsi_fusion/deuterium.bp",
      .mass = 3.344494690818129e-27,
      .charge = 1.0,
      .atomic_number = 1,
      .tracer = false
   },
   ParticleGroupSpec{
      .name = "neutrons",
      .filepath = "",
      .mass = 1.67492750056e-27,
      .charge = 0.0,
      .atomic_number = 0,
      .tracer = false
   },
   ParticleGroupSpec{
      .name = "helium3",
      .filepath = "",
      .mass = 5.008234522189299e-27,
      .charge = 2.0,
      .atomic_number = 2,
      .tracer = false
   },
   ParticleGroupSpec{
      .name = "tritium",
      .filepath = "",
      .mass = 5.008268858816166e-27,
      .charge = 1.0,
      .atomic_number = 1,
      .tracer = false
   },
   ParticleGroupSpec{
      .name = "protons",
      .filepath = "",
      .mass = 1.67262192595e-27,
      .charge = 1.0,
      .atomic_number = 1,
      .tracer = false
   },
   ParticleGroupSpec{
      .name = "electrons",
      .filepath = "/data/rlsi_fusion/electrons.bp",
      .mass = 9.1093837139e-31,
      .charge = -1.0,
      .atomic_number = 0,
      .tracer = false
   },
   ParticleGroupSpec{
      .name = "photons",
      .filepath = "",
      .mass = 0.0,
      .charge = 0.0,
      .atomic_number = 0,
      .tracer = false
   }
};

inline constexpr std::array<CollisionSpec, 1> collision_spec = {
   CollisionSpec{
      .group1 = "electrons",
      .group2 = "deuterium",
      .channels = {"radiation"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = false,
      .radiation = RadiationSpec{
         .product1 = "photons",
         .cross_section_file = "",
         .production_multiplier = 100000000.0,
         .min_energy = 1000.0,
         .max_energy = 10000000000.0,
         .reduce_electron_energy = true
      },

   }
};

/*---------------------------------------------------------------/
/-                      Metrics Parameters                      -/
/---------------------------------------------------------------*/
enum class MetricType { ParticleDump, ParticleDiag, ParticleEnergy, FieldDump, FieldEnergy };

inline constexpr auto metric_data_path = "/home/cepheid/TriForce/game_engine/data/rlsi_fusion";
inline constexpr std::array<MetricType, 1> metric_spec = {
	MetricType::ParticleDiag
};

#endif //PROGRAM_PARAM_HPP
