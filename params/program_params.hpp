#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

inline constexpr auto nThreads = 24;

inline constexpr auto Nx = 51zu;
inline constexpr auto Ny = 51zu;
inline constexpr auto Nz = 51zu;
inline constexpr auto NHalo = 2zu;

inline constexpr std::array x_range = {0.0, 0.005};
inline constexpr std::array y_range = {0.0, 0.005};
inline constexpr std::array z_range = {0.0, 0.005};

inline constexpr auto dx = 0.0001;
inline constexpr auto dy = 0.0001;
inline constexpr auto dz = 0.0001;

inline constexpr auto cfl   = 0.95;
inline constexpr auto dt    = 1.829541541469147e-13;
inline constexpr auto t_end = 7.318166165876588e-10;
inline constexpr auto Nt    = 100001zu;

inline constexpr auto save_interval = 1000zu;
inline constexpr auto interpolation_order = 2zu;

inline constexpr auto Ncx = 50zu;
inline constexpr auto Ncy = 50zu;
inline constexpr auto Ncz = 50zu;

inline constexpr auto PMLDepth    = 10zu;
inline constexpr auto PMLGrade    = 3.5;
inline constexpr auto PMLAlphaMax = 0.2;
//inline constexpr auto PMLKappaMax = 1.0;

#endif //PROGRAM_PARAM_HPP
