#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>
#include <cmath>

using compute_t = double;
constexpr compute_t operator""_fp(const long double x) { return static_cast<compute_t>(x); }

inline constexpr std::size_t nThreads = 32;

inline constexpr std::size_t Nx = 1501;
inline constexpr std::size_t Ny = 2;
inline constexpr std::size_t Nz = 1501;
inline constexpr std::size_t nHalo = 0zu;

inline constexpr std::array x_range = {-15.0e-6_fp, 15.0e-6_fp};
// inline constexpr std::array x_range = {-3.0e-6_fp, 12.0e-6_fp};
inline constexpr std::array y_range = {0.0_fp, 2.0e-8_fp}; // 40 nm
inline constexpr std::array z_range = {-15.0e-6_fp, 15.0e-6_fp};

inline constexpr auto dx = (x_range[1] - x_range[0]) / static_cast<compute_t>(Nx - 1zu);
inline constexpr auto dy = (y_range[1] - y_range[0]) / static_cast<compute_t>(Ny - 1zu);
inline constexpr auto dz = (z_range[1] - z_range[0]) / static_cast<compute_t>(Nz - 1zu);

inline constexpr auto Axy = dx * dy;
inline constexpr auto Axz = dx * dz;
inline constexpr auto Ayz = dy * dz;

inline constexpr auto cfl = 0.848_fp / 1.732050807568877_fp;
inline constexpr auto total_time = 3.0e-13_fp; // 300 fs
// inline constexpr auto total_time = 2.0e-14_fp; // 20 fs
inline constexpr auto dt = 4.0e-17_fp; // 0.04 fs
inline constexpr auto Nt = static_cast<int>(total_time / dt) + 1;

inline constexpr std::size_t save_interval = 100;

inline constexpr std::size_t Ncx = Nx - 1;
inline constexpr std::size_t Ncy = Ny - 1;
inline constexpr std::size_t Ncz = Nz - 1;

constexpr std::size_t get_cid(const std::size_t i, const std::size_t j, const std::size_t k) {
  return k + (Ncz * j) + (Ncz * Ncy * i);
}

#endif //PROGRAM_PARAM_HPP
