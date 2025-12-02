#ifndef EM_DATA_HPP
#define EM_DATA_HPP

#include "array_utils.hpp"
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
   vector_t psi;
   array_t b{};
   array_t c{};
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

   y0bc_t Exy0;
   y0bc_t Ezy0;
   y0bc_t Hxy0;
   y0bc_t Hzy0;

   y1bc_t Exy1;
   y1bc_t Ezy1;
   y1bc_t Hxy1;
   y1bc_t Hzy1;

   z0bc_t Exz0;
   z0bc_t Eyz0;
   z0bc_t Hxz0;
   z0bc_t Hyz0;

   z1bc_t Exz1;
   z1bc_t Eyz1;
   z1bc_t Hxz1;
   z1bc_t Hyz1;

   mdspan_t Exf{Ex.data(), {ex_ext, ex_stride}};
   mdspan_t Eyf{Ey.data(), {ey_ext, ey_stride}};
   mdspan_t Ezf{Ez.data(), {ez_ext, ez_stride}};

   mdspan_t Jxf{Jx.data(), {ex_ext, ex_stride}};
   mdspan_t Jyf{Jy.data(), {ey_ext, ey_stride}};
   mdspan_t Jzf{Jz.data(), {ez_ext, ez_stride}};
   
   mdspan_t Hxf{Hx.data(), {hx_ext, hx_stride}};
   mdspan_t Hyf{Hy.data(), {hy_ext, hy_stride}};
   mdspan_t Hzf{Hz.data(), {hz_ext, hz_stride}};

   mdspan_t Bxf{Hx.data(), {hx_ext, hx_stride}};
   mdspan_t Byf{Hy.data(), {hy_ext, hy_stride}};
   mdspan_t Bzf{Hz.data(), {hz_ext, hz_stride}};

   mdspan_t eyx0_psi{Eyx0.psi.data(), {eyhz_x_full_ext, eyhz_x_stride}};
   mdspan_t eyx1_psi{Eyx1.psi.data(), {eyhz_x_full_ext, eyhz_x_stride}};
   mdspan_t ezx0_psi{Ezx0.psi.data(), {hyez_x_full_ext, hyez_x_stride}};
   mdspan_t ezx1_psi{Ezx1.psi.data(), {hyez_x_full_ext, hyez_x_stride}};
   
   mdspan_t exy0_psi{Exy0.psi.data(), {exhz_y_full_ext, exhz_y_stride}};
   mdspan_t exy1_psi{Exy1.psi.data(), {exhz_y_full_ext, exhz_y_stride}};
   mdspan_t ezy0_psi{Ezy0.psi.data(), {hxez_y_full_ext, hxez_y_stride}};
   mdspan_t ezy1_psi{Ezy1.psi.data(), {hxez_y_full_ext, hxez_y_stride}};
   
   mdspan_t exz0_psi{Exz0.psi.data(), {exhy_z_full_ext, exhy_z_stride}};
   mdspan_t exz1_psi{Exz1.psi.data(), {exhy_z_full_ext, exhy_z_stride}};
   mdspan_t eyz0_psi{Eyz0.psi.data(), {hxey_z_full_ext, hxey_z_stride}};
   mdspan_t eyz1_psi{Eyz1.psi.data(), {hxey_z_full_ext, hxey_z_stride}};

   mdspan_t hyx0_psi{Hyx0.psi.data(), {hyez_x_full_ext, hyez_x_stride}};
   mdspan_t hyx1_psi{Hyx1.psi.data(), {hyez_x_full_ext, hyez_x_stride}};
   mdspan_t hzx0_psi{Hzx0.psi.data(), {eyhz_x_full_ext, eyhz_x_stride}};
   mdspan_t hzx1_psi{Hzx1.psi.data(), {eyhz_x_full_ext, eyhz_x_stride}};

   mdspan_t hxy0_psi{Hxy0.psi.data(), {hxez_y_full_ext, hxez_y_stride}};
   mdspan_t hxy1_psi{Hxy1.psi.data(), {hxez_y_full_ext, hxez_y_stride}};
   mdspan_t hzy0_psi{Hzy0.psi.data(), {exhz_y_full_ext, exhz_y_stride}};
   mdspan_t hzy1_psi{Hzy1.psi.data(), {exhz_y_full_ext, exhz_y_stride}};

   mdspan_t hxz0_psi{Hxz0.psi.data(), {hxey_z_full_ext, hxey_z_stride}};
   mdspan_t hxz1_psi{Hxz1.psi.data(), {hxey_z_full_ext, hxey_z_stride}};
   mdspan_t hyz0_psi{Hyz0.psi.data(), {exhy_z_full_ext, exhy_z_stride}};
   mdspan_t hyz1_psi{Hyz1.psi.data(), {exhy_z_full_ext, exhy_z_stride}};
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

   emdata.Exf = mdspan_t{emdata.Ex.data(), {ex_ext, ex_stride}};
   emdata.Eyf = mdspan_t{emdata.Ey.data(), {ey_ext, ey_stride}};
   emdata.Ezf = mdspan_t{emdata.Ez.data(), {ez_ext, ez_stride}};
   emdata.Jxf = mdspan_t{emdata.Jx.data(), {ex_ext, ex_stride}};
   emdata.Jyf = mdspan_t{emdata.Jy.data(), {ey_ext, ey_stride}};
   emdata.Jzf = mdspan_t{emdata.Jz.data(), {ez_ext, ez_stride}};
   emdata.Hxf = mdspan_t{emdata.Hx.data(), {hx_ext, hx_stride}};
   emdata.Hyf = mdspan_t{emdata.Hy.data(), {hy_ext, hy_stride}};
   emdata.Hzf = mdspan_t{emdata.Hz.data(), {hz_ext, hz_stride}};

   // constexpr auto e_coeff = dt / constants::eps0<double>;
   // constexpr auto h_coeff = dt / constants::mu0<double>;
   // constexpr auto e_max = std::max(ex_size, std::max(ey_size, ez_size));
   // constexpr auto h_max = std::max(hx_size, std::max(hy_size, hz_size));
   // emdata.Cexy = vector_t(e_max, e_coeff / dz);
   // emdata.Cexz = vector_t(e_max, e_coeff / dy);
   // emdata.Ceyz = vector_t(e_max, e_coeff / dx);
   // emdata.Chxy = vector_t(h_max, h_coeff / dz);
   // emdata.Chxz = vector_t(h_max, h_coeff / dy);
   // emdata.Chyz = vector_t(h_max, h_coeff / dx);
   // emdata.Chxy2 = vector_t(h_max, 0.5 * h_coeff / dz);
   // emdata.Chxz2 = vector_t(h_max, 0.5 * h_coeff / dy);
   // emdata.Chyz2 = vector_t(h_max, 0.5 * h_coeff / dx);

   init_bc<x0bc_t, X, Lo,  true>(emdata.Eyx0, std::extents{BCDepth, Ny - 1, Nz});
   init_bc<x0bc_t, X, Lo,  true>(emdata.Ezx0, std::extents{BCDepth, Ny, Nz - 1});

   init_bc<x0bc_t, X, Lo, false>(emdata.Hyx0, std::extents{BCDepth, Ny, Nz - 1});
   init_bc<x0bc_t, X, Lo, false>(emdata.Hzx0, std::extents{BCDepth, Ny - 1, Nz});

   init_bc<x1bc_t, X, Hi,  true>(emdata.Eyx1, std::extents{BCDepth, Ny - 1, Nz});
   init_bc<x1bc_t, X, Hi,  true>(emdata.Ezx1, std::extents{BCDepth, Ny, Nz - 1});

   init_bc<x1bc_t, X, Hi, false>(emdata.Hyx1, std::extents{BCDepth, Ny, Nz - 1});
   init_bc<x1bc_t, X, Hi, false>(emdata.Hzx1, std::extents{BCDepth, Ny - 1, Nz});

   init_bc<y0bc_t, Y, Lo,  true>(emdata.Exy0, std::extents{Nx - 1, BCDepth, Nz});
   init_bc<y0bc_t, Y, Lo,  true>(emdata.Ezy0, std::extents{Nx    , BCDepth, Nz - 1});

   init_bc<y0bc_t, Y, Lo, false>(emdata.Hxy0, std::extents{Nx    , BCDepth, Nz - 1});
   init_bc<y0bc_t, Y, Lo, false>(emdata.Hzy0, std::extents{Nx - 1, BCDepth, Nz});

   init_bc<y1bc_t, Y, Hi,  true>(emdata.Exy1, std::extents{Nx - 1, BCDepth, Nz});
   init_bc<y1bc_t, Y, Hi,  true>(emdata.Ezy1, std::extents{Nx    , BCDepth, Nz - 1});

   init_bc<y1bc_t, Y, Hi, false>(emdata.Hxy1, std::extents{Nx    , BCDepth, Nz - 1});
   init_bc<y1bc_t, Y, Hi, false>(emdata.Hzy1, std::extents{Nx - 1, BCDepth, Nz});

   init_bc<z0bc_t, Z, Lo,  true>(emdata.Exz0, std::extents{Nx - 1, Ny, BCDepth});
   init_bc<z0bc_t, Z, Lo,  true>(emdata.Eyz0, std::extents{Nx, Ny - 1, BCDepth});

   init_bc<z0bc_t, Z, Lo, false>(emdata.Hxz0, std::extents{Nx, Ny - 1, BCDepth});
   init_bc<z0bc_t, Z, Lo, false>(emdata.Hyz0, std::extents{Nx - 1, Ny, BCDepth});

   init_bc<z1bc_t, Z, Hi,  true>(emdata.Exz1, std::extents{Nx, Ny - 1, BCDepth});
   init_bc<z1bc_t, Z, Hi,  true>(emdata.Eyz1, std::extents{Nx - 1, Ny, BCDepth});

   init_bc<z1bc_t, Z, Hi, false>(emdata.Hxz1, std::extents{Nx, Ny - 1, BCDepth});
   init_bc<z1bc_t, Z, Hi, false>(emdata.Hyz1, std::extents{Nx - 1, Ny, BCDepth});

   emdata.eyx0_psi = mdspan_t{emdata.Eyx0.psi.data(), {eyhz_x_full_ext, eyhz_x_stride}};
   emdata.eyx1_psi = mdspan_t{emdata.Eyx1.psi.data(), {eyhz_x_full_ext, eyhz_x_stride}};
   emdata.ezx0_psi = mdspan_t{emdata.Ezx0.psi.data(), {hyez_x_full_ext, hyez_x_stride}};
   emdata.ezx1_psi = mdspan_t{emdata.Ezx1.psi.data(), {hyez_x_full_ext, hyez_x_stride}};
   
   emdata.exy0_psi = mdspan_t{emdata.Exy0.psi.data(), {exhz_y_full_ext, exhz_y_stride}};
   emdata.exy1_psi = mdspan_t{emdata.Exy1.psi.data(), {exhz_y_full_ext, exhz_y_stride}};
   emdata.ezy0_psi = mdspan_t{emdata.Ezy0.psi.data(), {hxez_y_full_ext, hxez_y_stride}};
   emdata.ezy1_psi = mdspan_t{emdata.Ezy1.psi.data(), {hxez_y_full_ext, hxez_y_stride}};
   
   emdata.exz0_psi = mdspan_t{emdata.Exz0.psi.data(), {exhy_z_full_ext, exhy_z_stride}};
   emdata.exz1_psi = mdspan_t{emdata.Exz1.psi.data(), {exhy_z_full_ext, exhy_z_stride}};
   emdata.eyz0_psi = mdspan_t{emdata.Eyz0.psi.data(), {hxey_z_full_ext, hxey_z_stride}};
   emdata.eyz1_psi = mdspan_t{emdata.Eyz1.psi.data(), {hxey_z_full_ext, hxey_z_stride}};

   emdata.hyx0_psi = mdspan_t{emdata.Hyx0.psi.data(), {hyez_x_full_ext, hyez_x_stride}};
   emdata.hyx1_psi = mdspan_t{emdata.Hyx1.psi.data(), {hyez_x_full_ext, hyez_x_stride}};
   emdata.hzx0_psi = mdspan_t{emdata.Hzx0.psi.data(), {eyhz_x_full_ext, eyhz_x_stride}};
   emdata.hzx1_psi = mdspan_t{emdata.Hzx1.psi.data(), {eyhz_x_full_ext, eyhz_x_stride}};
   
   emdata.hxy0_psi = mdspan_t{emdata.Hxy0.psi.data(), {hxez_y_full_ext, hxez_y_stride}};
   emdata.hxy1_psi = mdspan_t{emdata.Hxy1.psi.data(), {hxez_y_full_ext, hxez_y_stride}};
   emdata.hzy0_psi = mdspan_t{emdata.Hzy0.psi.data(), {exhz_y_full_ext, exhz_y_stride}};
   emdata.hzy1_psi = mdspan_t{emdata.Hzy1.psi.data(), {exhz_y_full_ext, exhz_y_stride}};
   
   emdata.hxz0_psi = mdspan_t{emdata.Hxz0.psi.data(), {hxey_z_full_ext, hxey_z_stride}};
   emdata.hxz1_psi = mdspan_t{emdata.Hxz1.psi.data(), {hxey_z_full_ext, hxey_z_stride}};
   emdata.hyz0_psi = mdspan_t{emdata.Hyz0.psi.data(), {exhy_z_full_ext, exhy_z_stride}};
   emdata.hyz1_psi = mdspan_t{emdata.Hyz1.psi.data(), {exhy_z_full_ext, exhy_z_stride}};
   
   return emdata;
} // end make_emdata()
} // end namespace tf::electromagnetics


#endif //EM_DATA_HPP
