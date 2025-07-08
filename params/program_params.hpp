#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include "compute_type.hpp"

#include <array>

inline constexpr auto nThreads = 48;

inline constexpr auto Nx = 1501zu;
inline constexpr auto Ny = 2zu;
inline constexpr auto Nz = 1501zu;

inline constexpr std::array x_range = {-1.5e-05_fp, 1.5e-05_fp};
inline constexpr std::array y_range = {0.0_fp, 2e-08_fp};
inline constexpr std::array z_range = {-1.5e-05_fp, 1.5e-05_fp};

inline constexpr auto dx = 2e-08_fp;
inline constexpr auto dy = 2e-08_fp;
inline constexpr auto dz = 2e-08_fp;

inline constexpr auto cfl   = 0.4497558596987185_fp;
inline constexpr auto dt    = 3e-17_fp;
inline constexpr auto t_end = 3e-13_fp;
inline constexpr auto Nt    = 10001zu;

inline constexpr auto save_interval = 100zu;
inline constexpr auto interpolation_order = 2zu;

inline constexpr auto Ncx = 1500zu;
inline constexpr auto Ncy = 1zu;
inline constexpr auto Ncz = 1500zu;

inline constexpr auto PMLDepth    = 20zu;
inline constexpr auto PMLGrade    = 3.5_fp;
inline constexpr auto PMLAlphaMax = 0.2_fp;
//inline constexpr auto PMLKappaMax = 1.0_fp;

#endif //PROGRAM_PARAM_HPP
