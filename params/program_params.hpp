#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>
#include <string>
#include <tuple>

inline constexpr auto nThreads = 24;

inline constexpr auto x_collapsed = true;
inline constexpr auto y_collapsed = true;
inline constexpr auto z_collapsed = true;

inline constexpr auto Nx = 2lu;
inline constexpr auto Ny = 2lu;
inline constexpr auto Nz = 2lu;

inline constexpr std::array x_range = {0.0, 1e-06};
inline constexpr std::array y_range = {0.0, 1e-06};
inline constexpr std::array z_range = {0.0, 1e-06};

inline constexpr auto dx = 1e-06;
inline constexpr auto dy = 1e-06;
inline constexpr auto dz = 1e-06;

inline constexpr auto cfl   = 0.2596278844909794;
inline constexpr auto dt    = 5e-16;
inline constexpr auto t_end = 5e-12;
inline constexpr auto Nt    = 10000lu;

inline constexpr auto sim_name = "carbon_thermal_eq";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = false;
inline constexpr auto push_enabled = false;
inline constexpr auto jdep_enabled = false;
inline constexpr auto coll_enabled = true;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
// enum class EMFace { X, Y, Z };
// enum class EMSide { Lo, Hi };

inline constexpr auto em_save_interval = 1lu;

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
enum class ParticleBCType { Static, Reflecting, Periodic, Outflow };

inline constexpr auto particle_save_interval = 100lu;
inline constexpr auto interpolation_order = 2lu;

inline constexpr auto PBCSelect = ParticleBCType::Outflow;

inline constexpr std::array particle_data = {"/data/carbon1.bp", "/data/carbon2.bp"};

inline constexpr std::array<collision_spec, 3> collision_params = {
   std::tuple("carbon1", "carbon2", 10, 1.0, 1, false),
   std::tuple("carbon1", "carbon1", 10, 1.0, 1, true),
   std::tuple("carbon2", "carbon2", 10, 1.0, 1, true),
};

/*---------------------------------------------------------------/
/-                      Metrics Parameters                      -/
/---------------------------------------------------------------*/
// inline constexpr std::array field_slices = {
//    {"Ex", 0, -1, -1} // Name, x, y, z (-1 is full extent)
// };

#endif //PROGRAM_PARAM_HPP
