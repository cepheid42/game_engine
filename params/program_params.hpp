#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

inline constexpr auto nThreads = 24;

inline constexpr auto Nx = 45zu;
inline constexpr auto Ny = 45zu;
inline constexpr auto Nz = 101zu;
inline constexpr auto NHalo = 2zu;

inline constexpr std::array x_range = {-0.11, 0.11};
inline constexpr std::array y_range = {-0.11, 0.11};
inline constexpr std::array z_range = {-0.25, 0.25};

inline constexpr auto dx = 0.005;
inline constexpr auto dy = 0.005;
inline constexpr auto dz = 0.005;

inline constexpr auto cfl   = 0.95;
inline constexpr auto dt    = 5e-12;
inline constexpr auto t_end = 8e-08;
inline constexpr auto Nt    = 16001zu;

inline constexpr auto save_interval = 40zu;
inline constexpr auto interpolation_order = 1zu;

inline constexpr auto Ncx = 44zu;
inline constexpr auto Ncy = 44zu;
inline constexpr auto Ncz = 100zu;

inline constexpr auto PMLDepth    = 6zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

#endif //PROGRAM_PARAM_HPP
