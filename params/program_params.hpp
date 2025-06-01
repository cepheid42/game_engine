#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "compute_type.hpp"

#include <array>

inline constexpr auto nThreads = 1;

inline constexpr auto Nx = 11zu;
inline constexpr auto Ny = 2zu;
inline constexpr auto Nz = 11zu;

inline constexpr std::array x_range = {0.0_fp, 8e-05_fp};
inline constexpr std::array y_range = {0.0_fp, 8e-06_fp};
inline constexpr std::array z_range = {0.0_fp, 8e-05_fp};

inline constexpr auto dx = 8.000000000000001e-06_fp;
inline constexpr auto dy = 8e-06_fp;
inline constexpr auto dz = 8.000000000000001e-06_fp;

inline constexpr auto cfl   = 0.5484827557301445_fp;
inline constexpr auto dt    = 1.4636e-14_fp;
inline constexpr auto t_end = 8e-10_fp;
inline constexpr auto Nt    = 54660zu;

inline constexpr auto save_interval = 1zu;

inline constexpr auto Ncx = Nx - 1zu;
inline constexpr auto Ncy = Ny - 1zu;
inline constexpr auto Ncz = Nz - 1zu;

inline constexpr auto PMLDepth    = 1zu;
inline constexpr auto PMLGrade    = 3.5_fp;
inline constexpr auto PMLAlphaMax = 0.2_fp;
//inline constexpr auto PMLKappaMax = 1.0_fp;

#endif //PROGRAM_PARAM_HPP
