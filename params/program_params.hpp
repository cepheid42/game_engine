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

inline constexpr auto dt    = 2e-17;
inline constexpr auto t_end = 3e-13;
inline constexpr auto Nt    = 15000zu;

inline constexpr auto sim_name = "lsi_full";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = true;
inline constexpr auto push_enabled = false;
inline constexpr auto jdep_enabled = false;
inline constexpr auto coll_enabled = false;
inline constexpr auto applied_fields_only = false;
inline constexpr auto velocity_backstep_enabled = false;
inline constexpr auto ionization_test_enabled = false;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 150zu;
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

constexpr auto laser_spec = LaserSpec{.lambda=8e-07, .E0=-27500000000000.0, .w0=2.5479e-06, .xspot=1.5e-05, .scale=1.28855495};

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Reflecting, Periodic, Outflow };
enum class ParticlePushType { Ballistic, Boris, HigueraCary };

inline constexpr auto particle_save_interval = 150zu;
inline constexpr auto sort_frequency = 100zu;
inline constexpr auto interpolation_order = 1zu;
inline constexpr auto ParticlePushSelect = ParticlePushType::Boris;
inline constexpr auto PBCSelect = ParticleBCType::Outflow;
inline constexpr auto PBCDepth = 3zu;

inline constexpr std::array<ParticleGroupSpec, 7> particle_spec = {
   ParticleGroupSpec{
      .name = "deuterium",
      .filepath = "/data/lsi_full/deuterium.bp",
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
      .filepath = "/data/lsi_full/electrons.bp",
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

inline constexpr std::array<CollisionSpec, 3> collision_spec = {
   CollisionSpec{
      .group1 = "electrons",
      .group2 = "electrons",
      .channels = {"coulomb"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = true,
      .coulomb = CoulombSpec{.coulomb_log = 0.0, .rate_multiplier = 1.0},
			
   },
   CollisionSpec{
      .group1 = "electrons",
      .group2 = "deuterium",
      .channels = {"coulomb", "radiation"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = false,
      .coulomb = CoulombSpec{.coulomb_log = 0.0, .rate_multiplier = 1.0},
			.radiation = RadiationSpec{
         .product1 = "photons",
         .cross_section_file = "/home/cepheid/TriForce/game_engine/tests/cross_section_data/SB_G4_Z1_kdsdk_MeV_barns.csv",
         .production_multiplier = 100000000.0,
         .min_energy = 0.0,
         .max_energy = 0.0,
         .reduce_electron_energy = true,
         .use_TFD = false
      },

   },
   CollisionSpec{
      .group1 = "deuterium",
      .group2 = "deuterium",
      .channels = {"coulomb", "fusion"},
      .step_interval = 1,
      .probability_search_area = 1.0,
      .self_scatter = true,
      .coulomb = CoulombSpec{.coulomb_log = 0.0, .rate_multiplier = 1.0},
		.fusion = {
         FusionSpec{
            .product1 = "neutrons",
            .product2 = "helium3",
            .cross_section_file = "/home/cepheid/TriForce/game_engine/tests/cross_section_data/DD_nHe3_BH_eV_m2.txt",
            .energy_gain = 3269000.0,
            .rate_multiplier = 1.0,
            .production_multiplier = 100000000.0,
            .constant_cross_section = 0.0
         },
         FusionSpec{
            .product1 = "tritium",
            .product2 = "protons",
            .cross_section_file = "/home/cepheid/TriForce/game_engine/tests/cross_section_data/DD_pT_BH_eV_m2.txt",
            .energy_gain = 4030000.0,
            .rate_multiplier = 1.0,
            .production_multiplier = 100000000.0,
            .constant_cross_section = 0.0
         },
      },
	
   }
};

/*---------------------------------------------------------------/
/-                      Metrics Parameters                      -/
/---------------------------------------------------------------*/
enum class MetricType { ParticleDump, ParticleDiag, ParticleEnergy, FieldDump, FieldEnergy };

inline constexpr auto metric_data_path = "/home/cepheid/TriForce/game_engine/data/lsi_full";
inline constexpr std::array<MetricType, 1> metric_spec = {
	MetricType::FieldDump
};

#endif //PROGRAM_PARAM_HPP
