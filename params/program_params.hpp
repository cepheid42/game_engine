#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "compute_type.hpp"

#include <array>

inline constexpr auto nThreads = 16;

inline constexpr auto Nx = 101zu;
inline constexpr auto Ny = 2zu;
inline constexpr auto Nz = 101zu;

inline constexpr std::array x_range = {-0.075_fp, 0.075_fp};
inline constexpr std::array y_range = {0.0_fp, 0.005_fp};
inline constexpr std::array z_range = {-0.075_fp, 0.075_fp};

inline constexpr auto dx = 0.0015_fp;
inline constexpr auto dy = 0.005_fp;
inline constexpr auto dz = 0.0015_fp;

inline constexpr auto cfl   = 0.14446830961840537_fp;
inline constexpr auto dt    = 5e-13_fp;
inline constexpr auto t_end = 5e-09_fp;
inline constexpr auto Nt    = 300001zu;

inline constexpr auto save_interval = 1000zu;
inline constexpr auto interpolation_order = 2zu;

inline constexpr auto Ncx = 100zu;
inline constexpr auto Ncy = 1zu;
inline constexpr auto Ncz = 100zu;

inline constexpr auto PMLDepth    = 1zu;
inline constexpr auto PMLGrade    = 3.5_fp;
inline constexpr auto PMLAlphaMax = 0.2_fp;
//inline constexpr auto PMLKappaMax = 1.0_fp;

#endif //PROGRAM_PARAM_HPP
