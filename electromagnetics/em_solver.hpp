#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

// #include "array_utils.hpp"
#include "constants.hpp"
// #include "bc_functors.hpp"
#include "em_params.hpp"
#include "em_sources.hpp"
#include "program_params.hpp"

#include <print>

namespace tf::electromagnetics {
template<Derivative D>
static auto diff(const auto& d, const auto i, const auto j, const auto k) {
   if      constexpr (D == Derivative::Dx) { return d[i + 1, j, k] - d[i, j, k]; }
   else if constexpr (D == Derivative::Dy) { return d[i, j + 1, k] - d[i, j, k]; }
   else if constexpr (D == Derivative::Dz) { return d[i, j, k + 1] - d[i, j, k]; }
}

template<Derivative D1, Derivative D2>
struct FieldIntegrator {
   static void apply(const mdspan_t& f, const mdspan_t& d1, const mdspan_t& d2, const auto& src,
                     const auto c_d1, const auto c_d2, const auto c_src)
   {
      for (std::size_t i = 0; i < f.extent(0); ++i) {
         for (std::size_t j = 0; j < f.extent(1); ++j) {
            for (std::size_t k = 0; k < f.extent(2); ++k) {
               const auto diff1   = c_d1 * diff<D1>(d1, i, j, k);
               const auto diff2   = c_d2 * diff<D2>(d2, i, j, k);
               const auto current = c_src * src[i, j, k];
               f[i, j, k] = f[i, j, k] + (diff1 - diff2) - current;
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
}; // end struct FieldIntegrator

template<Derivative D>
struct BCIntegrator {
   // static constexpr void apply(const auto&, const auto&, const auto&, const auto&) {}

   static void apply(const auto& psi, const auto& b, const auto& c, const auto& f, const auto& d, const auto cf)
   {
      std::size_t ipml;
      for (auto i = 0zu; i < psi.extent(0); ++i) {
         for (auto j = 0zu; j < psi.extent(1); ++j) {
            for (auto k = 0zu; k < psi.extent(2); ++k) {
               if      constexpr (D == Derivative::Dx) { ipml = i; }
               else if constexpr (D == Derivative::Dy) { ipml = j; }
               else if constexpr (D == Derivative::Dz) { ipml = k; }
               psi[i, j, k] = b[ipml] * psi[i, j, k] + c[ipml] * diff<D>(d, i, j, k);
               f[i, j, k] += cf * psi[i, j, k]; // cf comes with -1 baked in for += or -=
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
}; // end struct BCIntegrator

struct EMSolver {
   static constexpr auto dtdx_mu = dt / (constants::mu0<double> * dx);
   static constexpr auto dtdy_mu = dt / (constants::mu0<double> * dy);
   static constexpr auto dtdz_mu = dt / (constants::mu0<double> * dz);

   static constexpr auto dtdz_eps = dt / (constants::eps0<double> * dz);
   static constexpr auto dtdy_eps = dt / (constants::eps0<double> * dy);
   static constexpr auto dtdx_eps = dt / (constants::eps0<double> * dx);
   static constexpr auto   dt_eps = dt /  constants::eps0<double>;

   struct empty {
      static constexpr auto operator[](const auto, const auto, const auto) { return 0.0; }
   };
   
   static void advance(auto& emdata, const auto t, const auto n) requires(em_enabled) {
      auto ricker = [&n]()
      {
         constexpr auto Md = 2.0;
         constexpr auto ppw = 20.0;
         const auto alpha = math::SQR(constants::pi<double> * (cfl * static_cast<double>(n) / ppw - Md));
         return (1.0 - 2.0 * alpha) * std::exp(-alpha);
      };

      updateH(emdata);
      updateHBCs(emdata);

      // updateJBCs();
      // apply_srcs(t);

      emdata.Ezf[Nx / 2, Ny / 2, Nz / 2] += ricker();
      updateE(emdata);
      updateEBCs(emdata);

      // particle_correction(); // for the particles and shit
      // zero_currents();       // also for the particles, don't need last week's currents
   }

   static void advance(auto&, const auto) requires (!em_enabled) {}

   static void updateE(auto& emdata) {
      // const auto   Ex = mdspan_t{&emdata.Exf[0, 1, 1], {ex_update_ext, ex_stride}};
      // const auto   Jx = mdspan_t{&emdata.Jxf[0, 1, 1], {ex_update_ext, ex_stride}};
      // const auto Hzdy = mdspan_t{&emdata.Hzf[0, 0, 1], {ex_update_ext, hz_stride}};
      // const auto Hydz = mdspan_t{&emdata.Hyf[0, 1, 0], {ex_update_ext, hy_stride}};

      const auto   Ey = mdspan_t{&emdata.Eyf[1, 0, 1], {ey_update_ext, ey_stride}};
      const auto   Jy = mdspan_t{&emdata.Jyf[1, 0, 1], {ey_update_ext, ey_stride}};
      const auto Hxdz = mdspan_t{&emdata.Hxf[1, 0, 0], {ey_update_ext, hx_stride}};
      const auto Hzdx = mdspan_t{&emdata.Hzf[0, 0, 1], {ey_update_ext, hz_stride}};

      const auto   Ez = mdspan_t{&emdata.Ezf[1, 1, 0], {ez_update_ext, ez_stride}};
      const auto   Jz = mdspan_t{&emdata.Jzf[1, 1, 0], {ez_update_ext, ez_stride}};
      const auto Hydx = mdspan_t{&emdata.Hyf[0, 1, 0], {ez_update_ext, hy_stride}};
      const auto Hxdy = mdspan_t{&emdata.Hxf[1, 0, 0], {ez_update_ext, hx_stride}};

      FieldIntegrator<Derivative::Dy, Derivative::Dz>::apply(
         mdspan_t{&emdata.Exf[0, 1, 1], {ex_update_ext, ex_stride}},
         mdspan_t{&emdata.Hzf[0, 0, 1], {ex_update_ext, hz_stride}},
         mdspan_t{&emdata.Hyf[0, 1, 0], {ex_update_ext, hy_stride}},
         mdspan_t{&emdata.Jxf[0, 1, 1], {ex_update_ext, ex_stride}},
         dtdy_eps, dtdx_eps, dt_eps
      );

      FieldIntegrator<Derivative::Dz, Derivative::Dx>::apply(Ey, Hxdz, Hzdx, Jy, dtdz_eps, dtdx_eps, dt_eps);
      FieldIntegrator<Derivative::Dx, Derivative::Dy>::apply(Ez, Hydx, Hxdy, Jz, dtdx_eps, dtdy_eps, dt_eps);
   }

   static void updateH(auto& emdata) {
      FieldIntegrator<Derivative::Dz, Derivative::Dy>::apply(emdata.Hxf, emdata.Eyf, emdata.Ezf, empty{}, dtdz_mu, dtdy_mu, 0.0);
      FieldIntegrator<Derivative::Dx, Derivative::Dz>::apply(emdata.Hyf, emdata.Ezf, emdata.Exf, empty{}, dtdx_mu, dtdz_mu, 0.0);
      FieldIntegrator<Derivative::Dy, Derivative::Dx>::apply(emdata.Hzf, emdata.Exf, emdata.Eyf, empty{}, dtdy_mu, dtdx_mu, 0.0);
   }

   static void particle_correction(auto& emdata) {
      // std::ranges::copy(emdata.Hx, emdata.Bx.begin());
      // std::ranges::copy(emdata.Hy, emdata.By.begin());
      // std::ranges::copy(emdata.Hz, emdata.Bz.begin());
      //
      // FieldIntegrator<Derivative::Dz, Derivative::Dy>::apply(emdata.Bxf, emdata.Eyf, emdata.Ezf, empty{}, 0.5 * dtdz_mu, 0.5 * dtdy_mu, 0.0);
      // FieldIntegrator<Derivative::Dx, Derivative::Dz>::apply(emdata.Byf, emdata.Ezf, emdata.Exf, empty{}, 0.5 * dtdx_mu, 0.5 * dtdz_mu, 0.0);
      // FieldIntegrator<Derivative::Dy, Derivative::Dx>::apply(emdata.Bzf, emdata.Exf, emdata.Eyf, empty{}, 0.5 * dtdy_mu, 0.5 * dtdx_mu, 0.0);
   }

   static void updateEBCs(auto& emdata) {
      // Eyx0
      const auto eyx0 = mdspan_t{&emdata.Eyf[1, 0, 0], {eyhz_x_ext, ey_stride}};
      const auto hzx0 = mdspan_t{&emdata.Hzf[0, 0, 0], {eyhz_x_ext, hz_stride}};
      const auto eyx0_psi = mdspan_t{&emdata.eyx0_psi[1, 0, 0], {eyhz_x_ext, eyhz_x_stride}};
      const auto eyx0_b = std::span{emdata.Eyx0.b};
      const auto eyx0_c = std::span{emdata.Eyx0.c};
      BCIntegrator<Derivative::Dx>::apply(eyx0_psi, eyx0_b, eyx0_c, eyx0, hzx0, -dtdx_eps);

      // Eyx1
      const auto eyx1 = mdspan_t{&emdata.Eyf[Nx - BCDepth, 0, 0], {eyhz_x_ext, ey_stride}};
      const auto hzx1 = mdspan_t{&emdata.Hzf[Nx - BCDepth - 1, 0, 0], {eyhz_x_ext, hz_stride}};
      const auto eyx1_psi = mdspan_t{&emdata.eyx1_psi[0, 0, 0], {eyhz_x_ext, eyhz_x_stride}};
      const auto eyx1_b = std::span{emdata.Eyx1.b};
      const auto eyx1_c = std::span{emdata.Eyx1.c};
      BCIntegrator<Derivative::Dx>::apply(eyx1_psi, eyx1_b, eyx1_c, eyx1, hzx1, -dtdx_eps);

      // Ezx0
      const auto ezx0 = mdspan_t{&emdata.Ezf[1, 0, 0], {hyez_x_ext, ez_stride}};
      const auto hyx0 = mdspan_t{&emdata.Hyf[0, 0, 0], {hyez_x_ext, hy_stride}};
      const auto ezx0_psi = mdspan_t{&emdata.ezx0_psi[1, 0, 0], {hyez_x_ext, hyez_x_stride}};
      const auto ezx0_b = std::span{emdata.Ezx0.b};
      const auto ezx0_c = std::span{emdata.Ezx0.c};
      BCIntegrator<Derivative::Dx>::apply(ezx0_psi, ezx0_b, ezx0_c, ezx0, hyx0, dtdx_eps);

      // Ezx1
      const auto ezx1 = mdspan_t{&emdata.Ezf[Nx - BCDepth, 0, 0], {hyez_x_ext, ez_stride}};
      const auto hyx1 = mdspan_t{&emdata.Hyf[Nx - BCDepth - 1, 0, 0], {hyez_x_ext, hy_stride}};
      const auto ezx1_psi = mdspan_t{&emdata.ezx1_psi[0, 0, 0], {hyez_x_ext, hyez_x_stride}};
      const auto ezx1_b = std::span{emdata.Ezx1.b};
      const auto ezx1_c = std::span{emdata.Ezx1.c};
      BCIntegrator<Derivative::Dx>::apply(ezx1_psi, ezx1_b, ezx1_c, ezx1, hyx1, dtdx_eps);

      // Exy0
      const auto exy0 = mdspan_t{&emdata.Exf[0, 1, 0], {exhz_y_ext, ex_stride}};
      const auto hzy0 = mdspan_t{&emdata.Hzf[0, 0, 0], {exhz_y_ext, hz_stride}};
      const auto exy0_psi = mdspan_t{&emdata.exy0_psi[0, 1, 0], {exhz_y_ext, exhz_y_stride}};
      const auto exy0_b = std::span{emdata.Exy0.b};
      const auto exy0_c = std::span{emdata.Exy0.c};
      BCIntegrator<Derivative::Dy>::apply(exy0_psi, exy0_b, exy0_c, exy0, hzy0, dtdy_eps);

      // Exy1
      const auto exy1 = mdspan_t{&emdata.Exf[0, Ny - BCDepth, 0], {exhz_y_ext, ex_stride}};
      const auto hzy1 = mdspan_t{&emdata.Hzf[0, Ny - BCDepth - 1, 0], {exhz_y_ext, hz_stride}};
      const auto exy1_psi = mdspan_t{&emdata.exy1_psi[0, 0, 0], {exhz_y_ext, exhz_y_stride}};
      const auto exy1_b = std::span{emdata.Exy1.b};
      const auto exy1_c = std::span{emdata.Exy1.c};
      BCIntegrator<Derivative::Dy>::apply(exy1_psi, exy1_b, exy1_c, exy1, hzy1, dtdy_eps);

      // Ezy0
      const auto ezy0 = mdspan_t{&emdata.Ezf[0, 1, 0], {hxez_y_ext, ez_stride}};
      const auto hxy0 = mdspan_t{&emdata.Hxf[0, 0, 0], {hxez_y_ext, hx_stride}};
      const auto ezy0_psi = mdspan_t{&emdata.ezy0_psi[0, 1, 0], {hxez_y_ext, hxez_y_stride}};
      const auto ezy0_b = std::span{emdata.Ezy0.b};
      const auto ezy0_c = std::span{emdata.Ezy0.c};
      BCIntegrator<Derivative::Dy>::apply(ezy0_psi, ezy0_b, ezy0_c, ezy0, hxy0, -dtdy_eps);

      // Ezy1
      const auto ezy1 = mdspan_t{&emdata.Ezf[0, Ny - BCDepth, 0], {hxez_y_ext, ez_stride}};
      const auto hxy1 = mdspan_t{&emdata.Hxf[0, Ny - BCDepth - 1, 0], {hxez_y_ext, hx_stride}};
      const auto ezy1_psi = mdspan_t{&emdata.ezy1_psi[0, 0, 0], {hxez_y_ext, hxez_y_stride}};
      const auto ezy1_b = std::span{emdata.Ezy1.b};
      const auto ezy1_c = std::span{emdata.Ezy1.c};
      BCIntegrator<Derivative::Dy>::apply(ezy1_psi, ezy1_b, ezy1_c, ezy1, hxy1, -dtdy_eps);

      // Exz0
      const auto exz0 = mdspan_t{&emdata.Exf[0, 0, 1], {exhy_z_ext, ex_stride}};
      const auto hyz0 = mdspan_t{&emdata.Hyf[0, 0, 0], {exhy_z_ext, hy_stride}};
      const auto exz0_psi = mdspan_t{&emdata.exz0_psi[0, 0, 1], {exhy_z_ext, exhy_z_stride}};
      const auto exz0_b = std::span{emdata.Exz0.b};
      const auto exz0_c = std::span{emdata.Exz0.c};
      BCIntegrator<Derivative::Dz>::apply(exz0_psi, exz0_b, exz0_c, exz0, hyz0, -dtdz_eps);

      // Exz1
      const auto exz1 = mdspan_t{&emdata.Exf[0, 0, Nz - BCDepth], {exhy_z_ext, ex_stride}};
      const auto hyz1 = mdspan_t{&emdata.Hyf[0, 0, Nz - BCDepth - 1], {exhy_z_ext, hy_stride}};
      const auto exz1_psi = mdspan_t{&emdata.exz1_psi[0, 0, 0], {exhy_z_ext, exhy_z_stride}};
      const auto exz1_b = std::span{emdata.Exz1.b};
      const auto exz1_c = std::span{emdata.Exz1.c};
      BCIntegrator<Derivative::Dz>::apply(exz1_psi, exz1_b, exz1_c, exz1, hyz1, -dtdz_eps);

      // Eyz0
      const auto eyz0 = mdspan_t{&emdata.Eyf[0, 0, 1], {hxey_z_ext, ey_stride}};
      const auto hxz0 = mdspan_t{&emdata.Hxf[0, 0, 0], {hxey_z_ext, hx_stride}};
      const auto eyz0_psi = mdspan_t{&emdata.eyz0_psi[0, 0, 1], {hxey_z_ext, hxey_z_stride}};
      const auto eyz0_b = std::span{emdata.Eyz0.b};
      const auto eyz0_c = std::span{emdata.Eyz0.c};
      BCIntegrator<Derivative::Dz>::apply(eyz0_psi, eyz0_b, eyz0_c, eyz0, hxz0, dtdz_eps);

      // Eyz1
      const auto eyz1 = mdspan_t{&emdata.Eyf[0, 0, Nz - BCDepth], {hxey_z_ext, ey_stride}};
      const auto hxz1 = mdspan_t{&emdata.Hxf[0, 0, Nz - BCDepth - 1], {hxey_z_ext, hx_stride}};
      const auto eyz1_psi = mdspan_t{&emdata.eyz1_psi[0, 0, 0], {hxey_z_ext, hxey_z_stride}};
      const auto eyz1_b = std::span{emdata.Eyz1.b};
      const auto eyz1_c = std::span{emdata.Eyz1.c};
      BCIntegrator<Derivative::Dz>::apply(eyz1_psi, eyz1_b, eyz1_c, eyz1, hxz1, dtdz_eps);
   } // end updateEBCs()

   static void updateHBCs(auto& emdata) {
      // Hyx0
      const auto hyx0 = mdspan_t{&emdata.Hyf[0, 0, 0], {hyez_x_ext, hy_stride}};
      const auto ezx0 = mdspan_t{&emdata.Ezf[0, 0, 0], {hyez_x_ext, ez_stride}};
      const auto hyx0_psi = mdspan_t{&emdata.hyx0_psi[0, 0, 0], {hyez_x_ext, hyez_x_stride}};
      const auto hyx0_b = std::span{emdata.Hyx0.b};
      const auto hyx0_c = std::span{emdata.Hyx0.c};
      BCIntegrator<Derivative::Dx>::apply(hyx0_psi, hyx0_b, hyx0_c, hyx0, ezx0, dtdx_mu);

      // Hyx1
      const auto hyx1 = mdspan_t{&emdata.Hyf[Nx - BCDepth, 0, 0], {hyez_x_ext, hy_stride}};
      const auto ezx1 = mdspan_t{&emdata.Ezf[Nx - BCDepth, 0, 0], {hyez_x_ext, ez_stride}};
      const auto hyx1_psi = mdspan_t{&emdata.hyx1_psi[0, 0, 0], {hyez_x_ext, hyez_x_stride}};
      const auto hyx1_b = std::span{emdata.Hyx1.b};
      const auto hyx1_c = std::span{emdata.Hyx1.c};
      BCIntegrator<Derivative::Dx>::apply(hyx1_psi, hyx1_b, hyx1_c, hyx1, ezx1, dtdx_mu);

      // Hzx0
      const auto hzx0 = mdspan_t{&emdata.Hzf[0, 0, 0], {eyhz_x_ext, hz_stride}};
      const auto eyx0 = mdspan_t{&emdata.Eyf[0, 0, 0], {eyhz_x_ext, ey_stride}};
      const auto hzx0_psi = mdspan_t{&emdata.hzx0_psi[0, 0, 0], {eyhz_x_ext, eyhz_x_stride}};
      const auto hzx0_b = std::span{emdata.Hzx0.b};
      const auto hzx0_c = std::span{emdata.Hzx0.c};
      BCIntegrator<Derivative::Dx>::apply(hzx0_psi, hzx0_b, hzx0_c, hzx0, eyx0, -dtdx_mu);

      // Hzx1
      const auto hzx1 = mdspan_t{&emdata.Hzf[Nx - BCDepth, 0, 0], {eyhz_x_ext, hz_stride}};
      const auto eyx1 = mdspan_t{&emdata.Eyf[Nx - BCDepth, 0, 0], {eyhz_x_ext, ey_stride}};
      const auto hzx1_psi = mdspan_t{&emdata.hzx1_psi[0, 0, 0], {eyhz_x_ext, eyhz_x_stride}};
      const auto hzx1_b = std::span{emdata.Hzx1.b};
      const auto hzx1_c = std::span{emdata.Hzx1.c};
      BCIntegrator<Derivative::Dx>::apply(hzx1_psi, hzx1_b, hzx1_c, hzx1, eyx1, -dtdx_mu);

      // Hxy0
      const auto hxy0 = mdspan_t{&emdata.Hxf[0, 0, 0], {hxez_y_ext, hx_stride}};
      const auto ezy0 = mdspan_t{&emdata.Ezf[0, 0, 0], {hxez_y_ext, ez_stride}};
      const auto hxy0_psi = mdspan_t{&emdata.hxy0_psi[0, 0, 0], {hxez_y_ext, hxez_y_stride}};
      const auto hxy0_b = std::span{emdata.Hxy0.b};
      const auto hxy0_c = std::span{emdata.Hxy0.c};
      BCIntegrator<Derivative::Dy>::apply(hxy0_psi, hxy0_b, hxy0_c, hxy0, ezy0, -dtdy_mu);

      // Hxy1
      const auto hxy1 = mdspan_t{&emdata.Hxf[0, Ny - BCDepth, 0], {hxez_y_ext, hx_stride}};
      const auto ezy1 = mdspan_t{&emdata.Ezf[0, Ny - BCDepth, 0], {hxez_y_ext, ez_stride}};
      const auto hxy1_psi = mdspan_t{&emdata.hxy1_psi[0, 0, 0], {hxez_y_ext, hxez_y_stride}};
      const auto hxy1_b = std::span{emdata.Hxy1.b};
      const auto hxy1_c = std::span{emdata.Hxy1.c};
      BCIntegrator<Derivative::Dy>::apply(hxy1_psi, hxy1_b, hxy1_c, hxy1, ezy1, -dtdy_mu);

      // Hzy0
      const auto hzy0 = mdspan_t{&emdata.Hzf[0, 0, 0], {exhz_y_ext, hz_stride}};
      const auto exy0 = mdspan_t{&emdata.Exf[0, 0, 0], {exhz_y_ext, ex_stride}};
      const auto hzy0_psi = mdspan_t{&emdata.hzy0_psi[0, 0, 0], {exhz_y_ext, exhz_y_stride}};
      const auto hzy0_b = std::span{emdata.Hzy0.b};
      const auto hzy0_c = std::span{emdata.Hzy0.c};
      BCIntegrator<Derivative::Dy>::apply(hzy0_psi, hzy0_b, hzy0_c, hzy0, exy0, dtdy_mu);

      // Hzy1
      const auto hzy1 = mdspan_t{&emdata.Hzf[0, Ny - BCDepth, 0], {exhz_y_ext, hz_stride}};
      const auto exy1 = mdspan_t{&emdata.Exf[0, Ny - BCDepth, 0], {exhz_y_ext, ex_stride}};
      const auto hzy1_psi = mdspan_t{&emdata.hzy1_psi[0, 0, 0], {exhz_y_ext, exhz_y_stride}};
      const auto hzy1_b = std::span{emdata.Hzy1.b};
      const auto hzy1_c = std::span{emdata.Hzy1.c};
      BCIntegrator<Derivative::Dy>::apply(hzy1_psi, hzy1_b, hzy1_c, hzy1, exy1, dtdy_mu);

      // Hxz0
      const auto hxz0 = mdspan_t{&emdata.Hxf[0, 0, 0], {hxey_z_ext, hx_stride}};
      const auto eyz0 = mdspan_t{&emdata.Eyf[0, 0, 0], {hxey_z_ext, ey_stride}};
      const auto hxz0_psi = mdspan_t{&emdata.hxz0_psi[0, 0, 0], {hxey_z_ext, hxey_z_stride}};
      const auto hxz0_b = std::span{emdata.Hxz0.b};
      const auto hxz0_c = std::span{emdata.Hxz0.c};
      BCIntegrator<Derivative::Dz>::apply(hxz0_psi, hxz0_b, hxz0_c, hxz0, eyz0, dtdz_mu);

      // Hxz1
      const auto hxz1 = mdspan_t{&emdata.Hxf[0, 0, Nz - BCDepth], {hxey_z_ext, hx_stride}};
      const auto eyz1 = mdspan_t{&emdata.Eyf[0, 0, Nz - BCDepth], {hxey_z_ext, ey_stride}};
      const auto hxz1_psi = mdspan_t{&emdata.hxz1_psi[0, 0, 0], {hxey_z_ext, hxey_z_stride}};
      const auto hxz1_b = std::span{emdata.Hxz1.b};
      const auto hxz1_c = std::span{emdata.Hxz1.c};
      BCIntegrator<Derivative::Dz>::apply(hxz1_psi, hxz1_b, hxz1_c, hxz1, eyz1, dtdz_mu);

      // Hyz0
      const auto hyz0 = mdspan_t{&emdata.Hyf[0, 0, 0], {exhy_z_ext, hy_stride}};
      const auto exz0 = mdspan_t{&emdata.Exf[0, 0, 0], {exhy_z_ext, ex_stride}};
      const auto hyz0_psi = mdspan_t{&emdata.hyz0_psi[0, 0, 0], {exhy_z_ext, exhy_z_stride}};
      const auto hyz0_b = std::span{emdata.Hyz0.b};
      const auto hyz0_c = std::span{emdata.Hyz0.c};
      BCIntegrator<Derivative::Dz>::apply(hyz0_psi, hyz0_b, hyz0_c, hyz0, exz0, -dtdz_mu);

      // Hyz1
      const auto hyz1 = mdspan_t{&emdata.Hyf[0, 0, Nz - BCDepth], {exhy_z_ext, hy_stride}};
      const auto exz1 = mdspan_t{&emdata.Exf[0, 0, Nz - BCDepth], {exhy_z_ext, ex_stride}};
      const auto hyz1_psi = mdspan_t{&emdata.hyz1_psi[0, 0, 0], {exhy_z_ext, exhy_z_stride}};
      const auto hyz1_b = std::span{emdata.Hyz1.b};
      const auto hyz1_c = std::span{emdata.Hyz1.c};
      BCIntegrator<Derivative::Dz>::apply(hyz1_psi, hyz1_b, hyz1_c, hyz1, exz1, -dtdz_mu);
   } // end updateHBCs()

   // void updateJBCs() {
   //    // Only used for periodic BCs
   //    X0BC::Jy::apply(emdata.Jy, empty{}, empty{}, bcdata.x0.Jy);
   //    X0BC::Jz::apply(emdata.Jz, empty{}, empty{}, bcdata.x0.Jz);
   //
   //    Y0BC::Jx::apply(emdata.Jx, empty{}, empty{}, bcdata.y0.Jx);
   //    Y0BC::Jz::apply(emdata.Jz, empty{}, empty{}, bcdata.y0.Jz);
   //
   //    Z0BC::Jx::apply(emdata.Jx, empty{}, empty{}, bcdata.z0.Jx);
   //    Z0BC::Jy::apply(emdata.Jy, empty{}, empty{}, bcdata.z0.Jy);
   // }
   //
   //
   // void apply_srcs(const double t) {
   //    for (const auto& src: emdata.srcs) {
   //       src.apply(t);
   //    }
   //
   //    for (const auto& src: emdata.beams) {
   //       src.apply(t);
   //    }
   // }
   //
   // void zero_currents() {
   //    std::ranges::fill(emdata.Jx, 0.0);
   //    std::ranges::fill(emdata.Jy, 0.0);
   //    std::ranges::fill(emdata.Jz, 0.0);
   // }
}; // end struct EMSolver

} // end namespace tf::electromagnetics

#endif //EM_SOLVER_HPP
