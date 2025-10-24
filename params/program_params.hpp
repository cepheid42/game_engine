#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>
#include <string>
#include <tuple>

inline constexpr auto nThreads = 4;

inline constexpr auto is_2D_XZ = true;

inline constexpr auto Nx = 2lu;
inline constexpr auto Ny = 2lu;
inline constexpr auto Nz = 2lu;

inline constexpr std::array x_range = {0.0, 1e-06};
inline constexpr std::array y_range = {0.0, 1e-06};
inline constexpr std::array z_range = {0.0, 1e-06};

inline constexpr auto dx = 1e-06;
inline constexpr auto dy = 1e-06;
inline constexpr auto dz = 1e-06;

inline constexpr auto cfl   = 0.5192557689819588;
inline constexpr auto dt    = 1e-15;
inline constexpr auto t_end = 1e-11;
inline constexpr auto Nt    = 10000lu;

inline constexpr auto save_interval = 100lu;

inline constexpr auto sim_name = "carbon_thermal_eq";
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

inline constexpr auto PMLDepth    = 10lu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0lu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {2lu, 2lu, 2lu, 2lu, 2lu, 2lu};

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
using collision_spec = std::tuple<std::string, std::string, double, double, int, bool>;
enum class ParticleBCType { Periodic, Outflow };

inline constexpr auto interpolation_order = 2lu;

inline constexpr auto PBCSelect = ParticleBCType::Outflow;

inline constexpr std::array particle_data = {"/data/carbon1.bp", "/data/carbon2.bp"};

inline constexpr std::array<collision_spec, 3> collision_params = {
   std::tuple("carbon1", "carbon2", 10, 1.0, 1, false),
   std::tuple("carbon1", "carbon1", 10, 1.0, 1, true),
   std::tuple("carbon2", "carbon2", 10, 1.0, 1, true),
};

#endif //PROGRAM_PARAM_HPP
