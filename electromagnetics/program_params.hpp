#ifndef PROGRAM_PARAM_HPP
#define PROGRAM_PARAM_HPP

#include <array>

inline constexpr std::size_t Nx = 101zu;
inline constexpr std::size_t Ny = 101zu;
inline constexpr std::size_t Nz = 101zu;

inline constexpr std::array<double, 2> x_range = {0.0, 1.0};
inline constexpr std::array<double, 2> y_range = {0.0, 1.0};
inline constexpr std::array<double, 2> z_range = {0.0, 1.0};

inline constexpr double dx = (x_range[1] - x_range[0]) / static_cast<double>(Nx - 1zu);
inline constexpr double dy = (y_range[1] - y_range[0]) / static_cast<double>(Ny - 1zu);
inline constexpr double dz = (z_range[1] - z_range[0]) / static_cast<double>(Nz - 1zu);

inline constexpr double dt = 1.0e-12;
inline constexpr double cfl = 0.9;

#endif //PROGRAM_PARAM_HPP
