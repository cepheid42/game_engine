#ifndef EM_DATA_HPP
#define EM_DATA_HPP

#include "constants.hpp"
#include "em_params.hpp"
#include "math_utils.hpp"
#include "mdspan.hpp"
#include "program_params.hpp"
#include "traits.hpp"

#include <cmath>
#include <unordered_map>
#include <algorithm>

namespace tf::electromagnetics {

struct ReflectingData {};

struct PeriodicData {
   std::size_t numInterior;
   std::size_t hiIndex;
};

struct PMLData {
   using coeffs_t = std::array<double, PMLDepth>;
   vector_t psi;
   coeffs_t b{};
   coeffs_t c{};
};

struct emdata_t {
   vector_t Ex;
   vector_t Ey;
   vector_t Ez;
   vector_t Jx;
   vector_t Jy;
   vector_t Jz;
   vector_t Hx;
   vector_t Hy;
   vector_t Hz;
   vector_t Bx;
   vector_t By;
   vector_t Bz;
   
   // vector_t Cj;
   // vector_t Cexy;
   // vector_t Cexz;
   // vector_t Ceyz;
   // vector_t Chxy;
   // vector_t Chxz;
   // vector_t Chyz;
   // vector_t Chxy2;
   // vector_t Chxz2;
   // vector_t Chyz2;

   x0bc_t x0_Ey;
   x0bc_t x0_Ez;
   x0bc_t x0_Hy;
   x0bc_t x0_Hz;
   
   x1bc_t x1_Ey;
   x1bc_t x1_Ez;
   x1bc_t x1_Hy;
   x1bc_t x1_Hz;

   y0bc_t y0_Ex;
   y0bc_t y0_Ez;
   y0bc_t y0_Hx;
   y0bc_t y0_Hz;
   
   y1bc_t y1_Ex;
   y1bc_t y1_Ez;
   y1bc_t y1_Hx;
   y1bc_t y1_Hz;

   z0bc_t z0_Ex;
   z0bc_t z0_Ey;
   z0bc_t z0_Hx;
   z0bc_t z0_Hy;
   
   z1bc_t z1_Ex;
   z1bc_t z1_Ey;
   z1bc_t z1_Hx;
   z1bc_t z1_Hy;
}; // end struct emdata_t


template<typename T, EMFace, EMSide, bool>
requires std::same_as<T, ReflectingData>
auto init_bc(auto&, const auto&) {}

template<typename T, EMFace F, EMSide S, bool isE>
requires std::same_as<T, PMLData>
auto init_bc(auto& bc, const auto& dims) {
   constexpr auto delta = (F == EMFace::X ? dx : (F == EMFace::Y ? dy : dz));
   constexpr auto eta0      = constants::eta0<double>;
   constexpr auto pml_grade = static_cast<double>(PMLGrade);
   static constexpr auto m = 1.0; // alpha exponent
   static constexpr auto kappa_bc = 1.0;

   bc.psi = vector_t(dims.extent(0) * dims.extent(1) * dims.extent(2));
   auto d = math::linspace(1.0, 0.0, PMLDepth, false);

   if (!isE) {
      constexpr auto hstep = 1.0 / (2.0 * static_cast<double>(PMLDepth));
      for (auto& x: d) { x -= hstep; }
   }

   const auto sigma_max = 0.8 * (pml_grade + 1.0) / (delta * eta0);
   auto sigma(d);
   for (auto& x: sigma) { x = sigma_max * std::pow(x, PMLGrade); }

   auto alpha(d);
   for (auto& x: alpha) { x = PMLAlphaMax * std::pow(1.0 - x, m); }

   constexpr auto coef1 = -dt / constants::eps0<double>;
   for (auto i = 0lu; i < PMLDepth; i++) {
      bc.b[i] = std::exp(coef1 * ((sigma[i] / kappa_bc) + alpha[i]));
      bc.c[i] = (sigma[i] * (bc.b[i] - 1.0)) / (kappa_bc * (sigma[i] + (kappa_bc * alpha[i])));
   }

   if constexpr (is_HiSide<S>) {
      std::ranges::reverse(bc.b);
      std::ranges::reverse(bc.c);
   }
}

template<typename T, EMFace F, EMSide, bool>
requires std::same_as<T, PeriodicData>
auto init_bc(auto& bc, const auto& dims) {
   auto get_numInterior = [&] {
      if      constexpr (is_XFace<F>) { return dims.extent(0) - (2 * nHalo); }
      else if constexpr (is_YFace<F>) { return dims.extent(1) - (2 * nHalo); }
      else                            { return dims.extent(2) - (2 * nHalo); }
   };

   auto get_hiIndex = [&] {
      if      constexpr (is_XFace<F>) { return dims.extent(0) - 1 - nHalo; }
      else if constexpr (is_YFace<F>) { return dims.extent(1) - 1 - nHalo; }
      else                            { return dims.extent(2) - 1 - nHalo; }
   };

   bc.numInterior = get_numInterior();
   bc.hiIndex = get_hiIndex();
}

inline auto make_emdata() -> emdata_t {
   using enum EMFace;
   using enum EMSide;
   
   // constexpr auto e_coeff = dt / constants::eps0<double>;
   // constexpr auto h_coeff = dt / constants::mu0<double>;
   // constexpr auto e_max = std::max(ex_size, std::max(ey_size, ez_size));
   // constexpr auto h_max = std::max(hx_size, std::max(hy_size, hz_size));
   
   emdata_t emdata;
   emdata.Ex = vector_t(ex_size);
   emdata.Jx = vector_t(ex_size);

   emdata.Ey = vector_t(ey_size);
   emdata.Jy = vector_t(ey_size);

   emdata.Ez = vector_t(ez_size);
   emdata.Jz = vector_t(ez_size);

   emdata.Hx = vector_t(hx_size);
   emdata.Bx = vector_t(hx_size);

   emdata.Hy = vector_t(hy_size);
   emdata.By = vector_t(hy_size);

   emdata.Hz = vector_t(hz_size);
   emdata.Bz = vector_t(hz_size);

   // emdata.Cexy = vector_t(e_max, e_coeff / dz);
   // emdata.Cexz = vector_t(e_max, e_coeff / dy);
   // emdata.Ceyz = vector_t(e_max, e_coeff / dx);
   // emdata.Chxy = vector_t(h_max, h_coeff / dz);
   // emdata.Chxz = vector_t(h_max, h_coeff / dy);
   // emdata.Chyz = vector_t(h_max, h_coeff / dx);
   // emdata.Chxy2 = vector_t(h_max, 0.5 * h_coeff / dz);
   // emdata.Chxz2 = vector_t(h_max, 0.5 * h_coeff / dy);
   // emdata.Chyz2 = vector_t(h_max, 0.5 * h_coeff / dx);

   init_bc<x0bc_t, X, Lo,  true>(emdata.x0_Ey, eyx_ext{});
   init_bc<x0bc_t, X, Lo,  true>(emdata.x0_Ez, ezx_ext{});
   init_bc<x0bc_t, X, Lo, false>(emdata.x0_Hy, hyx_ext{});
   init_bc<x0bc_t, X, Lo, false>(emdata.x0_Hz, hzx_ext{});

   init_bc<x1bc_t, X, Hi,  true>(emdata.x1_Ey, eyx_ext{});
   init_bc<x1bc_t, X, Hi,  true>(emdata.x1_Ez, ezx_ext{});
   init_bc<x1bc_t, X, Hi, false>(emdata.x1_Hy, hyx_ext{});
   init_bc<x1bc_t, X, Hi, false>(emdata.x1_Hz, ezx_ext{});

   init_bc<y0bc_t, Y, Lo,  true>(emdata.y0_Ex, exy_ext{});
   init_bc<y0bc_t, Y, Lo,  true>(emdata.y0_Ez, ezy_ext{});
   init_bc<y0bc_t, Y, Lo, false>(emdata.y0_Hx, hxy_ext{});
   init_bc<y0bc_t, Y, Lo, false>(emdata.y0_Hz, hzy_ext{});

   init_bc<y1bc_t, Y, Hi,  true>(emdata.y1_Ex, exy_ext{});
   init_bc<y1bc_t, Y, Hi,  true>(emdata.y1_Ez, ezy_ext{});
   init_bc<y1bc_t, Y, Hi, false>(emdata.y1_Hx, hxy_ext{});
   init_bc<y1bc_t, Y, Hi, false>(emdata.y1_Hz, hzy_ext{});

   init_bc<z0bc_t, Z, Lo,  true>(emdata.z0_Ex, exz_ext{});
   init_bc<z0bc_t, Z, Lo,  true>(emdata.z0_Ey, eyz_ext{});
   init_bc<z0bc_t, Z, Lo, false>(emdata.z0_Hx, hxz_ext{});
   init_bc<z0bc_t, Z, Lo, false>(emdata.z0_Hy, hyz_ext{});

   init_bc<z1bc_t, Z, Hi,  true>(emdata.z1_Ex, exz_ext{});
   init_bc<z1bc_t, Z, Hi,  true>(emdata.z1_Ey, eyz_ext{});
   init_bc<z1bc_t, Z, Hi, false>(emdata.z1_Hx, hxz_ext{});
   init_bc<z1bc_t, Z, Hi, false>(emdata.z1_Hy, hyz_ext{});
   
   return emdata;
} // end make_emdata()

auto get_EMMap(const auto& emdata) -> std::unordered_map<std::string, mdspan_t> {
   return {{"Ex", {emdata.Ex.data(), std::extents{Nx - 1, Ny, Nz}}},
           {"Ey", {emdata.Ey.data(), std::extents{Nx, Ny - 1, Nz}}},
           {"Ez", {emdata.Ez.data(), std::extents{Nx, Ny, Nz - 1}}},
           {"Jx", {emdata.Jx.data(), std::extents{Nx - 1, Ny, Nz}}},
           {"Jy", {emdata.Jy.data(), std::extents{Nx, Ny - 1, Nz}}},
           {"Jz", {emdata.Jz.data(), std::extents{Nx, Ny, Nz - 1}}},
           {"Hx", {emdata.Hx.data(), std::extents{Nx, Ny - 1, Nz - 1}}},
           {"Hy", {emdata.Hy.data(), std::extents{Nx - 1, Ny, Nz - 1}}},
           {"Hz", {emdata.Hz.data(), std::extents{Nx - 1, Ny - 1, Nz}}}};
} // end get_EMMap()
} // end namespace tf::electromagnetics


#endif //EM_DATA_HPP
