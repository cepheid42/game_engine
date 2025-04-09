#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

using compute_t = float;

constexpr compute_t operator""_fp(const long double x) { return static_cast<compute_t>(x); }

inline constexpr std::size_t nThreads = 16;
//
// inline constexpr std::size_t Nx = 65;
// inline constexpr std::size_t Ny = 65;
// inline constexpr std::size_t Nz = 65;
//
// inline constexpr std::array x_range = {-0.1_fp, 0.1_fp};
// inline constexpr std::array y_range = {-0.1_fp, 0.1_fp};
// inline constexpr std::array z_range = {-0.1_fp, 0.1_fp};
//
// inline constexpr compute_t dx = (x_range[1] - x_range[0]) / static_cast<compute_t>(Nx - 1zu);
// inline constexpr compute_t dy = (y_range[1] - y_range[0]) / static_cast<compute_t>(Ny - 1zu);
// inline constexpr compute_t dz = (z_range[1] - z_range[0]) / static_cast<compute_t>(Nz - 1zu);
//
// inline constexpr compute_t Axy = dx * dy;
// inline constexpr compute_t Axz = dx * dz;
// inline constexpr compute_t Ayz = dy * dz;
//
// inline constexpr std::size_t nHalo = 2zu;
//
// inline constexpr auto cfl = 0.95_fp / 1.732050807568877_fp;
// inline constexpr auto Nt = 4000.0_fp;
// inline constexpr auto dt = cfl * dx / 299792458.0_fp;
// inline constexpr auto total_time = Nt * dt;
//
// inline constexpr std::size_t save_interval = 5;

inline constexpr std::size_t Nx = 16;
inline constexpr std::size_t Ny = 16;
inline constexpr std::size_t Nz = 16;
inline constexpr std::size_t nHalo = 2zu;

inline constexpr std::array x_range = {0.0_fp, 4.7e-6_fp};
inline constexpr std::array y_range = {0.0_fp, 4.7e-6_fp};
inline constexpr std::array z_range = {0.0_fp, 4.7e-6_fp};

inline constexpr compute_t dx = (x_range[1] - x_range[0]) / static_cast<compute_t>(Nx - 1zu);
inline constexpr compute_t dy = (y_range[1] - y_range[0]) / static_cast<compute_t>(Ny - 1zu);
inline constexpr compute_t dz = (z_range[1] - z_range[0]) / static_cast<compute_t>(Nz - 1zu);

inline constexpr compute_t Axy = dx * dy;
inline constexpr compute_t Axz = dx * dz;
inline constexpr compute_t Ayz = dy * dz;


inline constexpr auto cfl = 0.95_fp / 1.732050807568877_fp;
inline constexpr auto Nt = 10000.0_fp;
inline constexpr auto dt = cfl * dx / 299792458.0_fp;
inline constexpr auto total_time = Nt * dt;

inline constexpr std::size_t save_interval = 50;

inline constexpr std::size_t Ncx = Nx - 1;
inline constexpr std::size_t Ncy = Ny - 1;
inline constexpr std::size_t Ncz = Nz - 1;

constexpr std::size_t get_cid(const std::size_t i, const std::size_t j, const std::size_t k) {
  return k + (Ncz * j) + (Ncz * Ncy * i);
}

#endif //PROGRAM_PARAM_HPP
