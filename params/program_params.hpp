#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

using compute_t = float;

inline constexpr std::size_t nThreads = 32;

inline constexpr std::size_t Nx = 101zu;
inline constexpr std::size_t Ny = 101zu;
inline constexpr std::size_t Nz = 101zu;

inline constexpr std::array x_range = {0.0f, 1.0f};
inline constexpr std::array y_range = {0.0f, 1.0f};
inline constexpr std::array z_range = {0.0f, 1.0f};

inline constexpr compute_t dx = (x_range[1] - x_range[0]) / static_cast<compute_t>(Nx - 1zu);
inline constexpr compute_t dy = (y_range[1] - y_range[0]) / static_cast<compute_t>(Ny - 1zu);
inline constexpr compute_t dz = (z_range[1] - z_range[0]) / static_cast<compute_t>(Nz - 1zu);

inline constexpr auto cfl = 0.95f / 1.732050807568877f;
inline constexpr auto Nt = 4000.0f;
inline constexpr auto dt = cfl * dx / 299792458.0f;
inline constexpr auto total_time = Nt * dt;

inline constexpr std::size_t save_interval = 40;

#endif //PROGRAM_PARAM_HPP
