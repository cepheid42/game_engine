#ifndef GAME_ENGINE_EM_PARAMS_HPP
#define GAME_ENGINE_EM_PARAMS_HPP

#include "program_params.hpp"
#include "mdspan.hpp"
#include "traits.hpp"

#include <vector>

namespace tf::electromagnetics
{
// enum class EMFace { X, Y, Z };
// enum class EMSide { Lo, Hi };

static inline constexpr auto ex_size = (Nx - 1) * Ny * Nz;
static inline constexpr auto ey_size = Nx * (Ny - 1) * Nz;
static inline constexpr auto ez_size = Nx * Ny * (Nz - 1);
static inline constexpr auto hx_size = Nx * (Ny - 1) * (Nz - 1);
static inline constexpr auto hy_size = (Nx - 1) * Ny * (Nz - 1);
static inline constexpr auto hz_size = (Nx - 1) * (Ny - 1) * Nz;

using eyx_ext = std::extents<std::size_t, PMLDepth, Ny - 1, Nz>; // These get changed from python
using ezx_ext = std::extents<std::size_t, PMLDepth, Ny, Nz - 1>; // These get changed from python
using hyx_ext = ezx_ext;
using hzx_ext = eyx_ext;

using exy_ext = std::extents<std::size_t, Nx - 1, PMLDepth, Nz>; // These get changed from python
using ezy_ext = std::extents<std::size_t, Nx, PMLDepth, Nz - 1>; // These get changed from python
using hxy_ext = ezy_ext;
using hzy_ext = exy_ext;

using exz_ext = std::extents<std::size_t, Nx - 1, Ny, PMLDepth>; // These get changed from python
using eyz_ext = std::extents<std::size_t, Nx, Ny - 1, PMLDepth>; // These get changed from python
using hxz_ext = eyz_ext;
using hyz_ext = exz_ext;


static inline constexpr std::array ex_offsets = {0, 0, 1, 1, 1, 1};
static inline constexpr std::array ey_offsets = {1, 1, 0, 0, 1, 1};
static inline constexpr std::array ez_offsets = {1, 1, 1, 1, 0, 0};

using mdspan_t = std::mdspan<double, std::extents<std::size_t, 3>>;
using vector_t = std::vector<double>;

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