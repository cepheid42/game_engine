#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>
#include <string>
#include <tuple>

inline constexpr auto nThreads = 8;

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

inline constexpr auto cfl   = 0.0025962788449097936;
inline constexpr auto dt    = 5e-18;
inline constexpr auto t_end = 3.18e-15;
inline constexpr auto Nt    = 637zu;

inline constexpr auto sim_name = "ionization";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = false;
inline constexpr auto push_enabled = true;
inline constexpr auto jdep_enabled = false;
inline constexpr auto coll_enabled = true;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 1zu;

inline constexpr auto PMLDepth    = 0zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {2zu, 2zu, 2zu, 2zu, 2zu, 2zu};

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
using collision_spec = std::tuple<const char*, const char*, const char*, const char*, double, double, int, bool, double, const char*>;
enum class ParticleBCType { Static, Reflecting, Periodic, Outflow };

inline constexpr auto particle_save_interval = 1zu;
inline constexpr auto interpolation_order = 1zu;

inline constexpr auto PBCSelect = ParticleBCType::Periodic;

inline constexpr std::array particle_data = {"/data/electrons.bp", "/data/Al.bp", "/data/Al+.empty"};
inline constexpr auto particle_beam_file = "";
inline constexpr std::array<collision_spec, 1> collision_params = {
   std::tuple("electrons", "Al", "electrons", "Al+", 10, 1.0, 1, false, 5.9858, "/data/al0cs.txt"),
};

#endif //PROGRAM_PARAM_HPP
