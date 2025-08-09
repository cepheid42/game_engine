#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

inline constexpr auto nThreads = 8;

inline constexpr auto Nx = 301zu;
inline constexpr auto Ny = 2zu;
inline constexpr auto Nz = 301zu;
inline constexpr auto nHalo = 0zu;

inline constexpr std::array x_range = {-0.075, 0.075};
inline constexpr std::array y_range = {0.0, 0.0005};
inline constexpr std::array z_range = {-0.075, 0.075};

inline constexpr auto dx = 0.0005;
inline constexpr auto dy = 0.0005;
inline constexpr auto dz = 0.0005;

inline constexpr auto cfl   = 0.9346603841675257;
inline constexpr auto dt    = 9e-13;
inline constexpr auto t_end = 1e-07;
inline constexpr auto Nt    = 60000zu;

inline constexpr auto save_interval = 40zu;
inline constexpr auto interpolation_order = 2zu;

inline constexpr auto Ncx = 300zu;
inline constexpr auto Ncy = 1zu;
inline constexpr auto Ncz = 300zu;

inline constexpr auto PMLDepth    = 1zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

#endif //PROGRAM_PARAM_HPP
