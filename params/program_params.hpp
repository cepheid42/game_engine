#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>
#include <string>
#include <tuple>

inline constexpr auto nThreads = 24;

inline constexpr auto x_collapsed = false;
inline constexpr auto y_collapsed = false;
inline constexpr auto z_collapsed = false;

inline constexpr auto Nx = 101lu;
inline constexpr auto Ny = 101lu;
inline constexpr auto Nz = 101lu;

inline constexpr std::array x_range = {0.0, 1.0};
inline constexpr std::array y_range = {0.0, 1.0};
inline constexpr std::array z_range = {0.0, 1.0};

inline constexpr auto dx = 0.01;
inline constexpr auto dy = 0.01;
inline constexpr auto dz = 0.01;

inline constexpr auto cfl   = 0.95;
inline constexpr auto dt    = 1.829541541469147e-11;
inline constexpr auto t_end = 7.3181661658765885e-09;
inline constexpr auto Nt    = 401zu;

inline constexpr auto sim_name = "em_test";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

inline constexpr auto   em_enabled = true;
inline constexpr auto push_enabled = false;
inline constexpr auto jdep_enabled = false;
inline constexpr auto coll_enabled = false;

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class Derivative { Dx, Dy, Dz };

inline constexpr auto em_save_interval = 4lu;

inline constexpr auto PMLDepth    = 10lu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 0zu;

inline constexpr auto BCDepth = PMLDepth; // or nHalo

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {2zu, 2zu, 2zu, 2zu, 2zu, 2zu};

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
using collision_spec = std::tuple<std::string, std::string, double, double, int, bool>;
enum class ParticleBCType { Static, Reflecting, Periodic, Outflow };

inline constexpr auto particle_save_interval = 100zu;
inline constexpr auto interpolation_order = 2zu;

inline constexpr auto PBCSelect = ParticleBCType::Outflow;

inline constexpr std::array particle_data = {"/data/carbon1.bp", "/data/carbon2.bp"};
inline constexpr std::array<collision_spec, 3> collision_params = {
   std::tuple("carbon1", "carbon2", 10.0, 1.0, 1, false),
   std::tuple("carbon1", "carbon1", 10.0, 1.0, 1, true),
   std::tuple("carbon2", "carbon2", 10.0, 1.0, 1, true),
};

/*---------------------------------------------------------------/
/-                      Metrics Parameters                      -/
/---------------------------------------------------------------*/
// inline constexpr std::array field_slices = {
//    {"Ex", 0, -1, -1} // Name, x, y, z (-1 is full extent)
// };

#endif //PROGRAM_PARAM_HPP
