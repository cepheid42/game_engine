#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

using compute_t = float;

constexpr compute_t operator""_fp(const long double x) { return static_cast<compute_t>(x); }

inline constexpr std::size_t nThreads = 32;

inline constexpr std::size_t Nx = 65;
inline constexpr std::size_t Ny = 65;
inline constexpr std::size_t Nz = 65;

inline constexpr std::array x_range = {-0.1_fp, 0.1_fp};
inline constexpr std::array y_range = {-0.1_fp, 0.1_fp};
inline constexpr std::array z_range = {-0.1_fp, 0.1_fp};

inline constexpr compute_t dx = (x_range[1] - x_range[0]) / static_cast<compute_t>(Nx - 1zu);
inline constexpr compute_t dy = (y_range[1] - y_range[0]) / static_cast<compute_t>(Ny - 1zu);
inline constexpr compute_t dz = (z_range[1] - z_range[0]) / static_cast<compute_t>(Nz - 1zu);

inline constexpr compute_t Axy = dx * dy;
inline constexpr compute_t Axz = dx * dz;
inline constexpr compute_t Ayz = dy * dz;

inline constexpr std::size_t nHalo = 2zu;

inline constexpr auto cfl = 0.95_fp / 1.732050807568877_fp;
inline constexpr auto Nt = 1000.0_fp;
inline constexpr auto dt = cfl * dx / 299792458.0_fp;
inline constexpr auto total_time = Nt * dt;

inline constexpr std::size_t save_interval = 50;

#endif //PROGRAM_PARAM_HPP
