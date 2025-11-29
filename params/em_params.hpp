#ifndef GAME_ENGINE_EM_PARAMS_HPP
#define GAME_ENGINE_EM_PARAMS_HPP

#include "program_params.hpp"
#include "mdspan.hpp"

#include <vector>

namespace tf::electromagnetics
{
using vector_t = std::vector<double>;
using mdspan_t = std::mdspan<double, std::dextents<std::size_t, 3>, std::layout_stride>;

inline constexpr auto ex_size = (Nx - 1) * Ny * Nz;
inline constexpr auto ey_size = Nx * (Ny - 1) * Nz;
inline constexpr auto ez_size = Nx * Ny * (Nz - 1);
inline constexpr auto hx_size = Nx * (Ny - 1) * (Nz - 1);
inline constexpr auto hy_size = (Nx - 1) * Ny * (Nz - 1);
inline constexpr auto hz_size = (Nx - 1) * (Ny - 1) * Nz;

inline constexpr auto ex_ext = std::extents{Nx - 1, Ny, Nz};
inline constexpr auto ey_ext = std::extents{Nx, Ny - 1, Nz};
inline constexpr auto ez_ext = std::extents{Nx, Ny, Nz - 1};
inline constexpr auto hx_ext = std::extents{Nx, Ny - 1, Nz - 1};
inline constexpr auto hy_ext = std::extents{Nx - 1, Ny, Nz - 1};
inline constexpr auto hz_ext = std::extents{Nx - 1, Ny - 1, Nz};

inline constexpr auto ex_update_ext = std::extents{Nx - 1, Ny - 2, Nz - 2};
inline constexpr auto ey_update_ext = std::extents{Nx - 2, Ny - 1, Nz - 2};
inline constexpr auto ez_update_ext = std::extents{Nx - 2, Ny - 2, Nz - 1};

inline constexpr auto ex_stride = std::array{Ny * Nz, Nz, 1zu};
inline constexpr auto ey_stride = std::array{(Ny - 1) * Nz, Nz, 1zu};
inline constexpr auto ez_stride = std::array{Ny * (Nz - 1), Nz - 1, 1zu};
inline constexpr auto hx_stride = std::array{(Ny - 1) * (Nz - 1), Nz - 1, 1zu};
inline constexpr auto hy_stride = ez_stride;
inline constexpr auto hz_stride = ey_stride;

struct PMLData;
struct PeriodicData;
struct ReflectingData;

using x0bc_t = PMLData;
using x1bc_t = PMLData;
using y0bc_t = PMLData;
using y1bc_t = PMLData;
using z0bc_t = PMLData;
using z1bc_t = PMLData;

} // end namespace tf::electromagnetics

#endif //GAME_ENGINE_EM_PARAMS_HPP
