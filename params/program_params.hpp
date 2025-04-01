#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

using compute_t = float;

inline constexpr std::size_t nThreads = 32;


inline constexpr std::size_t Nx = 33;
inline constexpr std::size_t Ny = 41;
inline constexpr std::size_t Nz = 49;

inline constexpr std::array x_range = {-0.08f, 0.08f};
inline constexpr std::array y_range = {-0.1f, 0.1f};
inline constexpr std::array z_range = {-0.13f, 0.13f};

inline constexpr compute_t dx = (x_range[1] - x_range[0]) / static_cast<compute_t>(Nx - 1zu);
inline constexpr compute_t dy = (y_range[1] - y_range[0]) / static_cast<compute_t>(Ny - 1zu);
inline constexpr compute_t dz = (z_range[1] - z_range[0]) / static_cast<compute_t>(Nz - 1zu);

inline constexpr compute_t Axy = dx * dy;
inline constexpr compute_t Axz = dx * dz;
inline constexpr compute_t Ayz = dy * dz;

inline constexpr auto cfl = 0.5192557689819587f;// 0.95f / 1.732050807568877f;
inline constexpr auto Nt = 5000.0f;
inline constexpr auto dt = 5.0e-12f; //cfl * dx / 299792458.0f;
inline constexpr auto total_time = Nt * dt;

inline constexpr std::size_t save_interval = 5;

#endif //PROGRAM_PARAM_HPP
