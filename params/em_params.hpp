#ifndef GAME_ENGINE_EM_PARAMS_HPP
#define GAME_ENGINE_EM_PARAMS_HPP

#include "program_params.hpp"
#include "mdspan.hpp"

#include <vector>

namespace tf::electromagnetics
{
inline constexpr auto ex_size = (Nx - 1) * Ny * Nz;
inline constexpr auto ey_size = Nx * (Ny - 1) * Nz;
inline constexpr auto ez_size = Nx * Ny * (Nz - 1);
inline constexpr auto hx_size = Nx * (Ny - 1) * (Nz - 1);
inline constexpr auto hy_size = (Nx - 1) * Ny * (Nz - 1);
inline constexpr auto hz_size = (Nx - 1) * (Ny - 1) * Nz;

using vector_t = std::vector<double>;
using mdspan_t = std::mdspan<double, std::dextents<std::size_t, 3>, std::layout_stride>;

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
