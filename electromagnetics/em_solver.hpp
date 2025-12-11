#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "constants.hpp"
#include "em_data.hpp"
#include "em_params.hpp"
#include "program_params.hpp"

#include <cassert>
#include <print>

namespace tf::electromagnetics
{
template <Derivative D, bool back>
static auto diff(const auto& d, const auto i, const auto j, const auto k) {
   if constexpr (D == Derivative::Dx and !x_collapsed) {
      return back ? d[i, j, k] - d[i - 1, j, k] : d[i + 1, j, k] - d[i, j, k];
   }
   else if constexpr (D == Derivative::Dy and !y_collapsed) {
      return back ? d[i, j, k] - d[i, j - 1, k] : d[i, j + 1, k] - d[i, j, k];
   }
   else if constexpr (D == Derivative::Dz and !z_collapsed) {
      return back ? d[i, j, k] - d[i, j, k - 1] : d[i, j, k + 1] - d[i, j, k];
   }
   else { return 0.0; }
}

template <Derivative D1, Derivative D2, bool back, typename T=mdspan_t>
struct FieldUpdate {
   static void apply(const mdspan_t& f, const mdspan_t& d1, const mdspan_t& d2, const T& src,
                     const auto c_d1, const auto c_d2, const auto c_src) {
      #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = 0; i < f.extent(0); ++i) {
         for (std::size_t j = 0; j < f.extent(1); ++j) {
            for (std::size_t k = 0; k < f.extent(2); ++k) {
               const auto diff1 = c_d1 * diff<D1, back>(d1, i, j, k);
               const auto diff2 = c_d2 * diff<D2, back>(d2, i, j, k);
               const auto current = c_src * src[i, j, k];
               f[i, j, k] = f[i, j, k] + (diff1 - diff2) - current;
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
}; // end struct FieldUpdate

template <Derivative D, bool back>
struct BoundaryUpdate {
   static void apply(const auto&, const auto&, const auto&, const auto&, const auto&, const auto)
      requires (D == Derivative::Dy and y_collapsed) {
   }

   static void apply(const auto& f, const auto& d, const auto& psi, const auto& b, const auto& c, const auto cf) {
      #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (auto i = 0zu; i < psi.extent(0); ++i) {
         for (auto j = 0zu; j < psi.extent(1); ++j) {
            for (auto k = 0zu; k < psi.extent(2); ++k) {
               std::size_t ipml;
               if constexpr (D == Derivative::Dx) { ipml = i; }
               else if constexpr (D == Derivative::Dy) { ipml = j; }
               else if constexpr (D == Derivative::Dz) { ipml = k; }
               // else { assert(false); }
               psi[i, j, k] = b[ipml] * psi[i, j, k] + c[ipml] * diff<D, back>(d, i, j, k);
               f[i, j, k] += cf * psi[i, j, k]; // cf comes with -1 baked in for += or -=
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
}; // end struct BoundaryUpdate

struct EMSolver {
   using Derivative::Dx;
   using Derivative::Dy;
   using Derivative::Dz;

   static constexpr auto dtdx_mu = dt / (constants::mu0<double> * dx);
   static constexpr auto dtdy_mu = dt / (constants::mu0<double> * dy);
   static constexpr auto dtdz_mu = dt / (constants::mu0<double> * dz);

   static constexpr auto dtdx_eps = dt / (constants::eps0<double> * dx);
   static constexpr auto dtdy_eps = dt / (constants::eps0<double> * dy);
   static constexpr auto dtdz_eps = dt / (constants::eps0<double> * dz);
   static constexpr auto dt_eps = dt / constants::eps0<double>;

   struct empty {
      static constexpr auto operator[](const auto, const auto, const auto) { return 0.0; }
   };

   static void advance(auto&, const auto) requires (!em_enabled) {
   }

   static void advance(auto& emdata, const auto n) requires(em_enabled) {
      auto ricker = [&n]()
      {
         // constexpr auto Md = 2.0;
         constexpr auto ppw = 25.0;
         const auto alpha = math::SQR(constants::pi<double> * ((cfl * static_cast<double>(n) - 0.0) / ppw - 1.0));
         return (1.0 - 2.0 * alpha) * std::exp(-alpha);
      };

      updateH(emdata);
      updateHBCs(emdata);

      emdata.Ezf[Nx / 2, Ny / 2, Nz / 2] += 10.0 * ricker();
      // updateJBCs();
      // apply_srcs(emdata, static_cast<double>(n) * dt);

      updateE(emdata);
      updateEBCs(emdata);

      // particle_correction(emdata); // for the particles and shit
      // zero_currents(emdata);       // also for the particles, don't need last week's currents
   }

   static void updateE(auto& emdata) {
      FieldUpdate<Dy, Dz, true>::apply(
         {&emdata.Exf[0, 1, 1], {std::extents{Nx - 1, Ny - 2, Nz - 2}, ex_stride}},
         {&emdata.Hzf[0, 1, 1], {std::extents{Nx - 1, Ny - 2, Nz - 2}, hz_stride}},
         {&emdata.Hyf[0, 1, 1], {std::extents{Nx - 1, Ny - 2, Nz - 2}, hy_stride}},
         {&emdata.Jxf[0, 1, 1], {std::extents{Nx - 1, Ny - 2, Nz - 2}, ex_stride}},
         dtdy_eps, dtdx_eps
         , dt_eps
      );

      FieldUpdate<Dz, Dx, true>::apply(
         {&emdata.Eyf[1, 0, 1], {std::extents{Nx - 2, Ny - 1, Nz - 2}, ey_stride}},
         {&emdata.Hxf[1, 0, 1], {std::extents{Nx - 2, Ny - 1, Nz - 2}, hx_stride}},
         {&emdata.Hzf[1, 0, 1], {std::extents{Nx - 2, Ny - 1, Nz - 2}, hz_stride}},
         {&emdata.Jyf[1, 0, 1], {std::extents{Nx - 2, Ny - 1, Nz - 2}, ey_stride}},
         dtdz_eps, dtdx_eps
         , dt_eps
      );

      FieldUpdate<Dx, Dy, true>::apply(
         {&emdata.Ezf[1, 1, 0], {std::extents{Nx - 2, Ny - 2, Nz - 1}, ez_stride}},
         {&emdata.Hyf[1, 1, 0], {std::extents{Nx - 2, Ny - 2, Nz - 1}, hy_stride}},
         {&emdata.Hxf[1, 1, 0], {std::extents{Nx - 2, Ny - 2, Nz - 1}, hx_stride}},
         {&emdata.Jzf[1, 1, 0], {std::extents{Nx - 2, Ny - 2, Nz - 1}, ez_stride}},
         dtdx_eps, dtdy_eps
         , dt_eps
      );
   }

   static void updateH(auto& emdata) {
      FieldUpdate<Dz, Dy, false, empty>::apply(
         {&emdata.Hxf[0, 0, 0], {std::extents{Nx, Ny - 1, Nz - 1}, hx_stride}},
         {&emdata.Eyf[0, 0, 0], {std::extents{Nx, Ny - 1, Nz - 1}, ey_stride}},
         {&emdata.Ezf[0, 0, 0], {std::extents{Nx, Ny - 1, Nz - 1}, ez_stride}},
         empty{},
         dtdz_mu,
         dtdy_mu
         , 0.0
      );

      FieldUpdate<Dx, Dz, false, empty>::apply(
         {&emdata.Hyf[0, 0, 0], {std::extents{Nx - 1, Ny, Nz - 1}, hy_stride}},
         {&emdata.Ezf[0, 0, 0], {std::extents{Nx - 1, Ny, Nz - 1}, ez_stride}},
         {&emdata.Exf[0, 0, 0], {std::extents{Nx - 1, Ny, Nz - 1}, ex_stride}},
         empty{},
         dtdx_mu,
         dtdz_mu
         , 0.0
      );

      FieldUpdate<Dy, Dx, false, empty>::apply(
         {&emdata.Hzf[0, 0, 0], {std::extents{Nx - 1, Ny - 1, Nz}, hz_stride}},
         {&emdata.Exf[0, 0, 0], {std::extents{Nx - 1, Ny - 1, Nz}, ex_stride}},
         {&emdata.Eyf[0, 0, 0], {std::extents{Nx - 1, Ny - 1, Nz}, ey_stride}},
         empty{},
         dtdy_mu,
         dtdx_mu
         , 0.0
      );
   }

   static void updateEBCs(auto& emdata) {
      // Eyx0
      BoundaryUpdate<Dx, true>::apply(
         mdspan_t{&emdata.Eyf[1, 0, 0], {eyhz_x_ext, ey_stride}},
         mdspan_t{&emdata.Hzf[1, 0, 0], {eyhz_x_ext, hz_stride}},
         mdspan_t{&emdata.eyx0_psi[1, 0, 0], {eyhz_x_ext, eyhz_x_stride}},
         std::span{emdata.Eyx0.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Eyx0.c.begin() + 1, BCDepth - 1},
         -dtdx_eps
      );

      // Ezx0
      BoundaryUpdate<Dx, true>::apply(
         mdspan_t{&emdata.Ezf[1, 0, 0], {hyez_x_ext, ez_stride}},
         mdspan_t{&emdata.Hyf[1, 0, 0], {hyez_x_ext, hy_stride}},
         mdspan_t{&emdata.ezx0_psi[1, 0, 0], {hyez_x_ext, hyez_x_stride}},
         std::span{emdata.Ezx0.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Ezx0.c.begin() + 1, BCDepth - 1},
         dtdx_eps
      );

      // Eyx1
      BoundaryUpdate<Dx, true>::apply(
         mdspan_t{&emdata.Eyf[Nx - BCDepth, 0, 0], {eyhz_x_ext, ey_stride}},
         mdspan_t{&emdata.Hzf[Nx - BCDepth, 0, 0], {eyhz_x_ext, hz_stride}},
         mdspan_t{&emdata.eyx1_psi[0, 0, 0], {eyhz_x_ext, eyhz_x_stride}},
         std::span{emdata.Eyx1.b.begin(), BCDepth - 1},
         std::span{emdata.Eyx1.c.begin(), BCDepth - 1},
         -dtdx_eps
      );

      // Ezx1
      BoundaryUpdate<Dx, true>::apply(
         mdspan_t{&emdata.Ezf[Nx - BCDepth, 0, 0], {hyez_x_ext, ez_stride}},
         mdspan_t{&emdata.Hyf[Nx - BCDepth, 0, 0], {hyez_x_ext, hy_stride}},
         mdspan_t{&emdata.ezx1_psi[0, 0, 0], {hyez_x_ext, hyez_x_stride}},
         std::span{emdata.Ezx1.b.begin(), BCDepth - 1},
         std::span{emdata.Ezx1.c.begin(), BCDepth - 1},
         dtdx_eps
      );

      // Exy0
      BoundaryUpdate<Dy, true>::apply(
         mdspan_t{&emdata.Exf[0, 1, 0], {exhz_y_ext, ex_stride}},
         mdspan_t{&emdata.Hzf[0, 1, 0], {exhz_y_ext, hz_stride}},
         mdspan_t{&emdata.exy0_psi[0, 1, 0], {exhz_y_ext, exhz_y_stride}},
         std::span{emdata.Exy0.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Exy0.c.begin() + 1, BCDepth - 1},
         dtdy_eps
      );

      // Ezy0
      BoundaryUpdate<Dy, true>::apply(
         mdspan_t{&emdata.Ezf[0, 1, 0], {hxez_y_ext, ez_stride}},
         mdspan_t{&emdata.Hxf[0, 1, 0], {hxez_y_ext, hx_stride}},
         mdspan_t{&emdata.ezy0_psi[0, 1, 0], {hxez_y_ext, hxez_y_stride}},
         std::span{emdata.Ezy0.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Ezy0.c.begin() + 1, BCDepth - 1},
         -dtdy_eps
      );

      // Exy1
      BoundaryUpdate<Dy, true>::apply(
         mdspan_t{&emdata.Exf[0, Ny - BCDepth, 0], {exhz_y_ext, ex_stride}},
         mdspan_t{&emdata.Hzf[0, Ny - BCDepth, 0], {exhz_y_ext, hz_stride}},
         mdspan_t{&emdata.exy1_psi[0, 0, 0], {exhz_y_ext, exhz_y_stride}},
         std::span{emdata.Exy1.b.begin(), BCDepth - 1},
         std::span{emdata.Exy1.c.begin(), BCDepth - 1},
         dtdy_eps
      );

      // Ezy1
      BoundaryUpdate<Dy, true>::apply(
         mdspan_t{&emdata.Ezf[0, Ny - BCDepth, 0], {hxez_y_ext, ez_stride}},
         mdspan_t{&emdata.Hxf[0, Ny - BCDepth, 0], {hxez_y_ext, hx_stride}},
         mdspan_t{&emdata.ezy1_psi[0, 0, 0], {hxez_y_ext, hxez_y_stride}},
         std::span{emdata.Ezy1.b.begin(), BCDepth - 1},
         std::span{emdata.Ezy1.c.begin(), BCDepth - 1},
         -dtdy_eps
      );

      // Exz0
      BoundaryUpdate<Dz, true>::apply(
         mdspan_t{&emdata.Exf[0, 0, 1], {exhy_z_ext, ex_stride}},
         mdspan_t{&emdata.Hyf[0, 0, 1], {exhy_z_ext, hy_stride}},
         mdspan_t{&emdata.exz0_psi[0, 0, 1], {exhy_z_ext, exhy_z_stride}},
         std::span{emdata.Exz0.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Exz0.c.begin() + 1, BCDepth - 1},
         -dtdz_eps
      );

      // Eyz0
      BoundaryUpdate<Dz, true>::apply(
         mdspan_t{&emdata.Eyf[0, 0, 1], {hxey_z_ext, ey_stride}},
         mdspan_t{&emdata.Hxf[0, 0, 1], {hxey_z_ext, hx_stride}},
         mdspan_t{&emdata.eyz0_psi[0, 0, 1], {hxey_z_ext, hxey_z_stride}},
         std::span{emdata.Eyz0.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Eyz0.c.begin() + 1, BCDepth - 1},
         dtdz_eps
      );

      // Exz1
      BoundaryUpdate<Dz, true>::apply(
         mdspan_t{&emdata.Exf[0, 0, Nz - BCDepth], {exhy_z_ext, ex_stride}},
         mdspan_t{&emdata.Hyf[0, 0, Nz - BCDepth], {exhy_z_ext, hy_stride}},
         mdspan_t{&emdata.exz1_psi[0, 0, 0], {exhy_z_ext, exhy_z_stride}},
         std::span{emdata.Exz1.b.begin(), BCDepth - 1},
         std::span{emdata.Exz1.c.begin(), BCDepth - 1},
         -dtdz_eps
      );

      // Eyz1
      BoundaryUpdate<Dz, true>::apply(
         mdspan_t{&emdata.Eyf[0, 0, Nz - BCDepth], {hxey_z_ext, ey_stride}},
         mdspan_t{&emdata.Hxf[0, 0, Nz - BCDepth], {hxey_z_ext, hx_stride}},
         mdspan_t{&emdata.eyz1_psi[0, 0, 0], {hxey_z_ext, hxey_z_stride}},
         std::span{emdata.Eyz1.b.begin(), BCDepth - 1},
         std::span{emdata.Eyz1.c.begin(), BCDepth - 1},
         dtdz_eps
      );
   } // end updateEBCs()

   static void updateHBCs(auto& emdata) {
      // Hyx0
      BoundaryUpdate<Dx, false>::apply(
         mdspan_t{&emdata.Hyf[0, 0, 0], {hyez_x_ext, hy_stride}},
         mdspan_t{&emdata.Ezf[0, 0, 0], {hyez_x_ext, ez_stride}},
         mdspan_t{&emdata.hyx0_psi[0, 0, 0], {hyez_x_ext, hyez_x_stride}},
         std::span{emdata.Hyx0.b.begin(), BCDepth - 1},
         std::span{emdata.Hyx0.c.begin(), BCDepth - 1},
         dtdx_mu
      );

      // Hzx0
      BoundaryUpdate<Dx, false>::apply(
         mdspan_t{&emdata.Hzf[0, 0, 0], {eyhz_x_ext, hz_stride}},
         mdspan_t{&emdata.Eyf[0, 0, 0], {eyhz_x_ext, ey_stride}},
         mdspan_t{&emdata.hzx0_psi[0, 0, 0], {eyhz_x_ext, eyhz_x_stride}},
         std::span{emdata.Hzx0.b.begin(), BCDepth - 1},
         std::span{emdata.Hzx0.c.begin(), BCDepth - 1},
         -dtdx_mu
      );

      // Hyx1
      BoundaryUpdate<Dx, false>::apply(
         mdspan_t{&emdata.Hyf[Nx - BCDepth, 0, 0], {hyez_x_ext, hy_stride}},
         mdspan_t{&emdata.Ezf[Nx - BCDepth, 0, 0], {hyez_x_ext, ez_stride}},
         mdspan_t{&emdata.hyx1_psi[0, 0, 0], {hyez_x_ext, hyez_x_stride}},
         std::span{emdata.Hyx1.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Hyx1.c.begin() + 1, BCDepth - 1},
         dtdx_mu
      );

      // Hzx1
      BoundaryUpdate<Dx, false>::apply(
         mdspan_t{&emdata.Hzf[Nx - BCDepth, 0, 0], {eyhz_x_ext, hz_stride}},
         mdspan_t{&emdata.Eyf[Nx - BCDepth, 0, 0], {eyhz_x_ext, ey_stride}},
         mdspan_t{&emdata.hzx1_psi[0, 0, 0], {eyhz_x_ext, eyhz_x_stride}},
         std::span{emdata.Hzx1.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Hzx1.c.begin() + 1, BCDepth - 1},
         -dtdx_mu
      );

      // Hxy0
      BoundaryUpdate<Dy, false>::apply(
         mdspan_t{&emdata.Hxf[0, 0, 0], {hxez_y_ext, hx_stride}},
         mdspan_t{&emdata.Ezf[0, 0, 0], {hxez_y_ext, ez_stride}},
         mdspan_t{&emdata.hxy0_psi[0, 0, 0], {hxez_y_ext, hxez_y_stride}},
         std::span{emdata.Hxy0.b.begin(), BCDepth - 1},
         std::span{emdata.Hxy0.c.begin(), BCDepth - 1},
         -dtdy_mu
      );
      // Hzy0
      BoundaryUpdate<Dy, false>::apply(
         mdspan_t{&emdata.Hzf[0, 0, 0], {exhz_y_ext, hz_stride}},
         mdspan_t{&emdata.Exf[0, 0, 0], {exhz_y_ext, ex_stride}},
         mdspan_t{&emdata.hzy0_psi[0, 0, 0], {exhz_y_ext, exhz_y_stride}},
         std::span{emdata.Hzy0.b.begin(), BCDepth - 1},
         std::span{emdata.Hzy0.c.begin(), BCDepth - 1},
         dtdy_mu
      );

      // Hxy1
      BoundaryUpdate<Dy, false>::apply(
         mdspan_t{&emdata.Hxf[0, Ny - BCDepth, 0], {hxez_y_ext, hx_stride}},
         mdspan_t{&emdata.Ezf[0, Ny - BCDepth, 0], {hxez_y_ext, ez_stride}},
         mdspan_t{&emdata.hxy1_psi[0, 0, 0], {hxez_y_ext, hxez_y_stride}},
         std::span{emdata.Hxy1.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Hxy1.c.begin() + 1, BCDepth - 1},
         -dtdy_mu
      );

      // Hzy1
      BoundaryUpdate<Dy, false>::apply(
         mdspan_t{&emdata.Hzf[0, Ny - BCDepth, 0], {exhz_y_ext, hz_stride}},
         mdspan_t{&emdata.Exf[0, Ny - BCDepth, 0], {exhz_y_ext, ex_stride}},
         mdspan_t{&emdata.hzy1_psi[0, 0, 0], {exhz_y_ext, exhz_y_stride}},
         std::span{emdata.Hzy1.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Hzy1.c.begin() + 1, BCDepth - 1},
         dtdy_mu
      );

      // Hxz0
      BoundaryUpdate<Dz, false>::apply(
         mdspan_t{&emdata.Hxf[0, 0, 0], {hxey_z_ext, hx_stride}},
         mdspan_t{&emdata.Eyf[0, 0, 0], {hxey_z_ext, ey_stride}},
         mdspan_t{&emdata.hxz0_psi[0, 0, 0], {hxey_z_ext, hxey_z_stride}},
         std::span{emdata.Hxz0.b.begin(), BCDepth - 1},
         std::span{emdata.Hxz0.c.begin(), BCDepth - 1},
         dtdz_mu
      );

      // Hyz0
      BoundaryUpdate<Dz, false>::apply(
         mdspan_t{&emdata.Hyf[0, 0, 0], {exhy_z_ext, hy_stride}},
         mdspan_t{&emdata.Exf[0, 0, 0], {exhy_z_ext, ex_stride}},
         mdspan_t{&emdata.hyz0_psi[0, 0, 0], {exhy_z_ext, exhy_z_stride}},
         std::span{emdata.Hyz0.b.begin(), BCDepth - 1},
         std::span{emdata.Hyz0.c.begin(), BCDepth - 1},
         -dtdz_mu
      );

      // Hxz1
      BoundaryUpdate<Dz, false>::apply(
         mdspan_t{&emdata.Hxf[0, 0, Nz - BCDepth], {hxey_z_ext, hx_stride}},
         mdspan_t{&emdata.Eyf[0, 0, Nz - BCDepth], {hxey_z_ext, ey_stride}},
         mdspan_t{&emdata.hxz1_psi[0, 0, 0], {hxey_z_ext, hxey_z_stride}},
         std::span{emdata.Hxz1.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Hxz1.c.begin() + 1, BCDepth - 1},
         dtdz_mu
      );

      // Hyz1
      BoundaryUpdate<Dz, false>::apply(
         mdspan_t{&emdata.Hyf[0, 0, Nz - BCDepth], {exhy_z_ext, hy_stride}},
         mdspan_t{&emdata.Exf[0, 0, Nz - BCDepth], {exhy_z_ext, ex_stride}},
         mdspan_t{&emdata.hyz1_psi[0, 0, 0], {exhy_z_ext, exhy_z_stride}},
         std::span{emdata.Hyz1.b.begin() + 1, BCDepth - 1},
         std::span{emdata.Hyz1.c.begin() + 1, BCDepth - 1},
         -dtdz_mu
      );
   } // end updateHBCs()

   // static void particle_correction(auto& emdata) {
   //    std::ranges::copy(emdata.Hx, emdata.Bx.begin());
   //    std::ranges::copy(emdata.Hy, emdata.By.begin());
   //    std::ranges::copy(emdata.Hz, emdata.Bz.begin());
   //
   //    FieldUpdate<Dz, Dy, empty>::apply(emdata.Bxf, emdata.Eyf, emdata.Ezf, empty{}, 0.5 * dtdz_mu, 0.5 * dtdy_mu, 0.0);
   //    FieldUpdate<Dx, Dz, empty>::apply(emdata.Byf, emdata.Ezf, emdata.Exf, empty{}, 0.5 * dtdx_mu, 0.5 * dtdz_mu, 0.0);
   //    FieldUpdate<Dy, Dx, empty>::apply(emdata.Bzf, emdata.Exf, emdata.Eyf, empty{}, 0.5 * dtdy_mu, 0.5 * dtdx_mu, 0.0);
   //
   //    std::ranges::for_each(emdata.Bx, [](auto& b){ b *= constants::mu0<double>; });
   //    std::ranges::for_each(emdata.By, [](auto& b){ b *= constants::mu0<double>; });
   //    std::ranges::for_each(emdata.Bz, [](auto& b){ b *= constants::mu0<double>; });
   //
   //    for (std::size_t i = 0; i < emdata.Bx.size(); i++) {
   //       emdata.Bx[i] = emdata.Bx[i] * constants::mu0<double>;
   //    }
   //
   //    for (std::size_t i = 0; i < emdata.By.size(); i++) {
   //       emdata.By[i] = emdata.By[i] * constants::mu0<double>;
   //    }
   //
   //    for (std::size_t i = 0; i < emdata.Bz.size(); i++) {
   //       emdata.Bz[i] = emdata.Bz[i] * constants::mu0<double>;
   //    }
   //
   //    // todo: add up total fields here
   // }

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


   // static void apply_srcs(auto& emdata, const auto t) {
   //    // for (const auto& src: emdata.srcs) {
   //    //    src.apply(t);
   //    // }
   //
   //    for (const auto& src: emdata.beams) {
   //       src.apply(t);
   //    }
   // }
   //
   // static void zero_currents(auto& emdata) {
   //    std::ranges::fill(emdata.Jx, 0.0);
   //    std::ranges::fill(emdata.Jy, 0.0);
   //    std::ranges::fill(emdata.Jz, 0.0);
   // }
}; // end struct EMSolver
} // end namespace tf::electromagnetics

#endif //EM_SOLVER_HPP
