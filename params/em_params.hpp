#ifndef GAME_ENGINE_EM_PARAMS_HPP
#define GAME_ENGINE_EM_PARAMS_HPP

#include "program_params.hpp"
#include "mdspan.hpp"

#include <array>
#include <span>
#include <vector>

namespace tf::electromagnetics
{
using vector_t = std::vector<double>;
using array_t = std::array<double, PMLDepth>;
using span_t = std::span<double, BCDepth>;
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

inline constexpr auto eyhz_x_full_ext = std::extents{BCDepth, Ny - 1, Nz    };
inline constexpr auto hyez_x_full_ext = std::extents{BCDepth, Ny    , Nz - 1};

inline constexpr auto exhz_y_full_ext = std::extents{Nx - 1, BCDepth, Nz    };
inline constexpr auto hxez_y_full_ext = std::extents{Nx    , BCDepth, Nz - 1};

inline constexpr auto exhy_z_full_ext = std::extents{Nx - 1, Ny    , BCDepth};
inline constexpr auto hxey_z_full_ext = std::extents{Nx    , Ny - 1, BCDepth};

inline constexpr auto eyhz_x_ext = std::extents{BCDepth - 1, Ny - 1, Nz    };
inline constexpr auto hyez_x_ext = std::extents{BCDepth - 1, Ny    , Nz - 1};

inline constexpr auto exhz_y_ext = std::extents{Nx - 1, BCDepth - 1, Nz    };
inline constexpr auto hxez_y_ext = std::extents{Nx    , BCDepth - 1, Nz - 1};

inline constexpr auto exhy_z_ext = std::extents{Nx - 1, Ny    , BCDepth - 1};
inline constexpr auto hxey_z_ext = std::extents{Nx    , Ny - 1, BCDepth - 1};

inline constexpr auto eyhz_x_stride = ey_stride;
inline constexpr auto hyez_x_stride = ez_stride;
inline constexpr auto exhz_y_stride = std::array{ BCDepth *       Nz,      Nz, 1zu};
inline constexpr auto hxez_y_stride = std::array{ BCDepth * (Nz - 1),  Nz - 1, 1zu};
inline constexpr auto exhy_z_stride = std::array{     Ny  *  BCDepth, BCDepth, 1zu};
inline constexpr auto hxey_z_stride = std::array{(Ny - 1) *  BCDepth, BCDepth, 1zu};

} // end namespace tf::electromagnetics

#endif //GAME_ENGINE_EM_PARAMS_HPP
