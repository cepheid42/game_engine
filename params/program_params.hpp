#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

inline constexpr auto nThreads = 32;

inline constexpr auto is_2D_XZ = true;

inline constexpr auto Nx = 1501lu;
inline constexpr auto Ny = 2lu;
inline constexpr auto Nz = 1501lu;

inline constexpr std::array x_range = {-1.5e-05, 1.5e-05};
inline constexpr std::array y_range = {0.0, 2e-08};
inline constexpr std::array z_range = {-1.5e-05, 1.5e-05};

inline constexpr auto dx = 2e-08;
inline constexpr auto dy = 2e-08;
inline constexpr auto dz = 2e-08;

inline constexpr auto cfl   = 1.0385115379639174;
inline constexpr auto dt    = 4e-17;
inline constexpr auto t_end = 3e-13;
inline constexpr auto Nt    = 7500lu;

inline constexpr auto save_interval = 100lu;

inline constexpr auto sim_name = "lsi";
inline constexpr auto sim_path = "/home/cepheid/TriForce/game_engine";

/*---------------------------------------------------------------/
/-                        EM Parameters                         -/
/---------------------------------------------------------------*/
enum class EMFace { X, Y, Z };
enum class EMSide { Lo, Hi };

inline constexpr auto PMLDepth    = 15lu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

inline constexpr auto nHalo = 2lu;

// Periodic = 0, PML = 1, Reflecting = 2
inline constexpr std::array BCSelect = {1lu, 1lu, 2lu, 2lu, 1lu, 1lu};

/*---------------------------------------------------------------/
/-                     Particle Parameters                      -/
/---------------------------------------------------------------*/
enum class ParticleBCType { Periodic, Outflow };

inline constexpr auto interpolation_order = 2lu;

inline constexpr auto PBCSelect = ParticleBCType::Outflow;

inline constexpr std::array particle_data = {"/data/electrons.bp", "/data/ions.bp"};
#endif //PROGRAM_PARAM_HPP
