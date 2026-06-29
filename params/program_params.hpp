#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "particle_spec.hpp"
#include "sources_spec.hpp"

#include <array>

inline constexpr auto nThreads = 4;

inline constexpr auto x_collapsed = true;
inline constexpr auto y_collapsed = true;
inline constexpr auto z_collapsed = true;

inline constexpr auto Nx = 2zu;
inline constexpr auto Ny = 2zu;
inline constexpr auto Nz = 2zu;

inline constexpr std::array x_range = {0.0, 1.6e-06};
inline constexpr std::array y_range = {0.0, 1e-06};
inline constexpr std::array z_range = {0.0, 1.6e-06};

inline constexpr auto dx = 1.6e-06;
inline constexpr auto dy = 1e-06;
inline constexpr auto dz = 1.6e-06;

inline constexpr auto dt    = 1e-15;
inline constexpr auto t_end = 4e-14;
inline constexpr auto Nt    = 41zu;

inline constexpr auto sim_name = "boron_brems";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = false;
inline constexpr auto push_enabled = false;
inline constexpr auto jdep_enabled = false;
inline constexpr auto coll_enabled = true;
inline constexpr auto applied_fields_only = false;
inline constexpr auto velocity_backstep_enabled = false;
inline constexpr auto ionization_test_enabled = false;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 10zu;
inline constexpr auto em_subcycles = 1zu;
inline constexpr auto dt_em = 1e-15;

inline constexpr auto PMLDepth    = 10zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {2zu, 2zu, 2zu, 2zu, 2zu, 2zu};

inline constexpr auto laser_enabled = false;
inline constexpr auto applied_fields_path = "";

constexpr auto laser_spec = LaserSpec{.lambda=0.0, .E0=0.0, .w0=0.0, .xspot=0.0, .scale=0.0};

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Reflecting, Periodic, Outflow };
enum class ParticlePushType { Ballistic, Boris, HigueraCary };

inline constexpr auto particle_save_interval = 1zu;
inline constexpr auto sort_frequency = 100zu;
inline constexpr auto interpolation_order = 1zu;
inline constexpr auto ParticlePushSelect = ParticlePushType::Boris;
inline constexpr auto PBCSelect = ParticleBCType::Periodic;
inline constexpr auto PBCDepth = 3zu;

inline constexpr std::array<ParticleGroupSpec, 3> particle_spec = {
   ParticleGroupSpec{
      .name = "electrons",
      .filepath = "/data/boron_brems/electrons.bp",
      .mass = 9.1093837139e-31,
      .charge = -1.0,
      .atomic_number = 0,
      .tracer = false
   },
   ParticleGroupSpec{
      .name = "boron",
      .filepath = "/data/boron_brems/boron.bp",
      .mass = 1.7952087874094118e-26,
      .charge = 1.0,
      .atomic_number = 5,
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
      .group2 = "boron",
      .channels = {"radiation"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = false,
      .radiation = RadiationSpec{
         .product1 = "photons",
         .cross_section_file = "/home/cepheid/TriForce/game_engine/tests/cross_section_data/SB_G4_Z5_kdsdk_MeV_barns.csv",
         .production_multiplier = 1e+30,
         .reduce_electron_energy = false
      },
	
   }
};

/*---------------------------------------------------------------/
/-                      Metrics Parameters                      -/
/---------------------------------------------------------------*/
enum class MetricType { ParticleDump, ParticleDiag, ParticleEnergy, FieldDump, FieldEnergy };

inline constexpr auto metric_data_path = "/home/cepheid/TriForce/game_engine/data/boron_brems";
inline constexpr std::array<MetricType, 1> metric_spec = {
	MetricType::ParticleDump
};

#endif //PROGRAM_PARAM_HPP
