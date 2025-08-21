#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

inline constexpr auto nThreads = 32;

inline constexpr auto Nx = 101zu;
inline constexpr auto Ny = 101zu;
inline constexpr auto Nz = 101zu;
inline constexpr auto NHalo = 0zu;

inline constexpr std::array x_range = {-0.01, 0.01};
inline constexpr std::array y_range = {-0.01, 0.01};
inline constexpr std::array z_range = {-0.01, 0.01};

inline constexpr auto dx = 0.0002;
inline constexpr auto dy = 0.0002;
inline constexpr auto dz = 0.0002;

inline constexpr auto cfl   = 0.95;
inline constexpr auto dt    = 3.659083082938294e-13;
inline constexpr auto t_end = 1.4636332331753175e-10;
inline constexpr auto Nt    = 400zu;

inline constexpr auto save_interval = 4zu;
inline constexpr auto interpolation_order = 2zu;

inline constexpr auto Ncx = 100zu;
inline constexpr auto Ncy = 100zu;
inline constexpr auto Ncz = 100zu;

inline constexpr auto PMLDepth    = 10zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

#endif //PROGRAM_PARAM_HPP
