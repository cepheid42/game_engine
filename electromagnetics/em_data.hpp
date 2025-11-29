#ifndef EM_DATA_HPP
#define EM_DATA_HPP

#include "constants.hpp"
#include "em_params.hpp"
#include "math_utils.hpp"
#include "mdspan.hpp"
#include "program_params.hpp"
#include "traits.hpp"

#include <cmath>
// #include <unordered_map>
#include <algorithm>

namespace tf::electromagnetics
{

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

   x0bc_t Eyx0;
   x0bc_t Ezx0;
   x0bc_t Hyx0;
   x0bc_t Hzx0;
   
   x1bc_t Eyx1;
   x1bc_t Ezx1;
   x1bc_t Hyx1;
   x1bc_t Hzx1;

   // y0bc_t Exy0;
   // y0bc_t Ezy0;
   // y0bc_t Hxy0;
   // y0bc_t Hzy0;
   //
   // y1bc_t Exy1;
   // y1bc_t Ezy1;
   // y1bc_t Hxy1;
   // y1bc_t Hzy1;
   //
   // z0bc_t Exz0;
   // z0bc_t Eyz0;
   // z0bc_t Hxz0;
   // z0bc_t Hyz0;
   //
   // z1bc_t Exz1;
   // z1bc_t Eyz1;
   // z1bc_t Hxz1;
   // z1bc_t Hyz1;
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

// template<typename T, EMFace F, EMSide, bool>
// requires std::same_as<T, PeriodicData>
// auto init_bc(auto& bc, const auto& dims) {
//    auto get_numInterior = [&] {
//       if      constexpr (is_XFace<F>) { return dims.extent(0) - (2 * nHalo); }
//       else if constexpr (is_YFace<F>) { return dims.extent(1) - (2 * nHalo); }
//       else                            { return dims.extent(2) - (2 * nHalo); }
//    };
//
//    auto get_hiIndex = [&] {
//       if      constexpr (is_XFace<F>) { return dims.extent(0) - 1 - nHalo; }
//       else if constexpr (is_YFace<F>) { return dims.extent(1) - 1 - nHalo; }
//       else                            { return dims.extent(2) - 1 - nHalo; }
//    };
//
//    bc.numInterior = get_numInterior();
//    bc.hiIndex = get_hiIndex();
// }

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

   init_bc<x0bc_t, X, Lo,  true>(emdata.Eyx0, eyx0_ext);
   init_bc<x0bc_t, X, Lo,  true>(emdata.Ezx0, ezx0_ext);
   // init_bc<x0bc_t, X, Lo, false>(emdata.Hyx0, hyx_ext);
   // init_bc<x0bc_t, X, Lo, false>(emdata.Hzx0, hzx_ext);

   init_bc<x1bc_t, X, Hi,  true>(emdata.Eyx1, eyx0_ext);
   init_bc<x1bc_t, X, Hi,  true>(emdata.Ezx1, ezx0_ext);
   // init_bc<x1bc_t, X, Hi, false>(emdata.Hyx1, hyx_ext{});
   // init_bc<x1bc_t, X, Hi, false>(emdata.Hzx1, ezx_ext{});

   // init_bc<y0bc_t, Y, Lo,  true>(emdata.Exy0, exy_ext{});
   // init_bc<y0bc_t, Y, Lo,  true>(emdata.Ezy0, ezy_ext{});
   // init_bc<y0bc_t, Y, Lo, false>(emdata.Hxy0, hxy_ext{});
   // init_bc<y0bc_t, Y, Lo, false>(emdata.Hzy0, hzy_ext{});
   //
   // init_bc<y1bc_t, Y, Hi,  true>(emdata.Exy1, exy_ext{});
   // init_bc<y1bc_t, Y, Hi,  true>(emdata.Ezy1, ezy_ext{});
   // init_bc<y1bc_t, Y, Hi, false>(emdata.Hxy1, hxy_ext{});
   // init_bc<y1bc_t, Y, Hi, false>(emdata.Hzy1, hzy_ext{});
   //
   // init_bc<z0bc_t, Z, Lo,  true>(emdata.Exz0, exz_ext{});
   // init_bc<z0bc_t, Z, Lo,  true>(emdata.Eyz0, eyz_ext{});
   // init_bc<z0bc_t, Z, Lo, false>(emdata.Hxz0, hxz_ext{});
   // init_bc<z0bc_t, Z, Lo, false>(emdata.Hyz0, hyz_ext{});
   //
   // init_bc<z1bc_t, Z, Hi,  true>(emdata.Exz1, exz_ext{});
   // init_bc<z1bc_t, Z, Hi,  true>(emdata.Eyz1, eyz_ext{});
   // init_bc<z1bc_t, Z, Hi, false>(emdata.Hxz1, hxz_ext{});
   // init_bc<z1bc_t, Z, Hi, false>(emdata.Hyz1, hyz_ext{});
   
   return emdata;
} // end make_emdata()
} // end namespace tf::electromagnetics


#endif //EM_DATA_HPP
