#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "constants.hpp"
// #include "bc_functors.hpp"
#include "em_params.hpp"
#include "em_sources.hpp"

#include <print>
#include <cassert>

namespace tf::electromagnetics {

template<Derivative D1, Derivative D2>
struct FieldIntegrator {
   template<Derivative D>
   static auto diff(const auto& d, const auto i, const auto j, const auto k) {
      if      constexpr (D == Derivative::FDx) { return d[i + 1, j, k] - d[i, j, k]; }
      else if constexpr (D == Derivative::FDy) { return d[i, j + 1, k] - d[i, j, k]; }
      else if constexpr (D == Derivative::FDz) { return d[i, j, k + 1] - d[i, j, k]; }
      // else if constexpr (D == Derivative::BDx) { return d[i, j, k] - d[i - 1, j, k]; }
      // else if constexpr (D == Derivative::BDy) { return d[i, j, k] - d[i, j - 1, k]; }
      // else if constexpr (D == Derivative::BDz) { return d[i, j, k] - d[i, j, k - 1]; }
   }

   static void apply(auto& f, const auto& d1, const auto& d2, const auto& src,
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


struct EMSolver {
   struct empty {
      static auto operator[](const auto, const auto, const auto) { return 0.0; }
   };
   
   static void advance(auto& emdata, const auto t, const auto n) requires(em_enabled) {
      auto ricker = [&n]()
      {
         constexpr auto Md = 2.0;
         constexpr auto ppw = 20.0;
         const auto alpha = math::SQR(constants::pi<double> * (cfl * static_cast<double>(n) / ppw - Md));
         return (1.0 - 2.0 * alpha) * std::exp(-alpha);
      };

      // const auto Ex = mdspan_t{emdata.Ex.data(), Nx - 1, Ny, Nz};
      // const auto Ey = mdspan_t{emdata.Ey.data(), Nx, Ny - 1, Nz};
      const auto Ez = mdspan_t{emdata.Ez.data(), {std::extents{Nx, Ny, Nz - 1}, std::array{Ny * (Nz - 1), Nz - 1, 1zu}}};

      updateH(emdata);
      // updateHBCs(emdata);
      // updateJBCs();
      // apply_srcs(t);
      Ez[Nx / 2, Ny / 2, Nz / 2] += ricker();
      updateE(emdata);
      // updateEBCs(emdata);
      // particle_correction(); // for the particles and shit
      // zero_currents();       // also for the particles, don't need last week's currents
   }

   static void advance(auto&, const auto) requires (!em_enabled) {}

   static void updateE(auto& emdata) {
      const auto Exf = mdspan_t{emdata.Ex.data(), {std::extents{Nx - 1, Ny, Nz}, std::array{Ny * Nz, Nz, 1zu}}};
      const auto Eyf = mdspan_t{emdata.Ey.data(), {std::extents{Nx, Ny - 1, Nz}, std::array{(Ny - 1) * Nz, Nz, 1zu}}};
      const auto Ezf = mdspan_t{emdata.Ez.data(), {std::extents{Nx, Ny, Nz - 1}, std::array{Ny * (Nz - 1), Nz - 1, 1zu}}};
      const auto Jxf = mdspan_t{emdata.Jx.data(), {std::extents{Nx - 1, Ny, Nz}, std::array{Ny * Nz, Nz, 1zu}}};
      const auto Jyf = mdspan_t{emdata.Jy.data(), {std::extents{Nx, Ny - 1, Nz}, std::array{(Ny - 1) * Nz, Nz, 1zu}}};
      const auto Jzf = mdspan_t{emdata.Jz.data(), {std::extents{Nx, Ny, Nz - 1}, std::array{Ny * (Nz - 1), Nz - 1, 1zu}}};
      const auto Hxf = mdspan_t{emdata.Hx.data(), {std::extents{Nx, Ny - 1, Nz - 1}, std::array{(Ny - 1) * (Nz - 1), Nz - 1, 1zu}}};
      const auto Hyf = mdspan_t{emdata.Hy.data(), {std::extents{Nx - 1, Ny, Nz - 1}, std::array{Ny * (Nz - 1), Nz - 1, 1zu}}};
      const auto Hzf = mdspan_t{emdata.Hz.data(), {std::extents{Nx - 1, Ny - 1, Nz}, std::array{(Ny - 1) * Nz, Nz, 1zu}}};

      const auto   Ex = mdspan_t{&Exf[0, 1, 1], {std::extents{Nx - 1, Ny - 2, Nz - 2}, std::array{ Ny      *  Nz     , Nz    , 1zu}}};
      const auto   Jx = mdspan_t{&Jxf[0, 1, 1], {std::extents{Nx - 1, Ny - 2, Nz - 2}, std::array{ Ny      *  Nz     , Nz    , 1zu}}};
      const auto Hzdy = mdspan_t{&Hzf[0, 0, 1], {std::extents{Nx - 1, Ny - 2, Nz - 2}, std::array{(Ny - 1) *  Nz     , Nz    , 1zu}}};
      const auto Hydz = mdspan_t{&Hyf[0, 1, 0], {std::extents{Nx - 1, Ny - 2, Nz - 2}, std::array{ Ny      * (Nz - 1), Nz - 1, 1zu}}};

      const auto   Ey = mdspan_t{&Eyf[1, 0, 1], {std::extents{Nx - 2, Ny - 1, Nz - 2}, std::array{(Ny - 1) *  Nz     , Nz    , 1zu}}};
      const auto   Jy = mdspan_t{&Jyf[1, 0, 1], {std::extents{Nx - 2, Ny - 1, Nz - 2}, std::array{(Ny - 1) *  Nz     , Nz    , 1zu}}};
      const auto Hxdz = mdspan_t{&Hxf[1, 0, 0], {std::extents{Nx - 2, Ny - 1, Nz - 2}, std::array{(Ny - 1) * (Nz - 1), Nz - 1, 1zu}}};
      const auto Hzdx = mdspan_t{&Hzf[0, 0, 1], {std::extents{Nx - 2, Ny - 1, Nz - 2}, std::array{(Ny - 1) *  Nz     , Nz    , 1zu}}};

      const auto   Ez = mdspan_t{&Ezf[1, 1, 0], {std::extents{Nx - 2, Ny - 2, Nz - 1}, std::array{ Ny      * (Nz - 1), Nz - 1, 1zu}}};
      const auto   Jz = mdspan_t{&Jzf[1, 1, 0], {std::extents{Nx - 2, Ny - 2, Nz - 1}, std::array{ Ny      * (Nz - 1), Nz - 1, 1zu}}};
      const auto Hydx = mdspan_t{&Hxf[0, 1, 0], {std::extents{Nx - 2, Ny - 2, Nz - 2}, std::array{(Ny - 1) * (Nz - 1), Nz - 1, 1zu}}};
      const auto Hxdy = mdspan_t{&Hzf[1, 0, 0], {std::extents{Nx - 2, Ny - 2, Nz - 1}, std::array{(Ny - 1) *  Nz     , Nz    , 1zu}}};

      constexpr auto Cexy = dt / (constants::eps0<double> * dz);
      constexpr auto Cexz = dt / (constants::eps0<double> * dy);
      constexpr auto Ceyz = dt / (constants::eps0<double> * dx);
      constexpr auto Cj = dt / constants::eps0<double>;

      FieldIntegrator<Derivative::FDy, Derivative::FDz>::apply(Ex, Hzdy, Hydz, Jx, Cexz, Ceyz, Cj);
      FieldIntegrator<Derivative::FDz, Derivative::FDx>::apply(Ey, Hxdz, Hzdx, Jy, Cexy, Ceyz, Cj);
      FieldIntegrator<Derivative::FDx, Derivative::FDy>::apply(Ez, Hydx, Hxdy, Jz, Ceyz, Cexz, Cj);
   }

   static void updateH(auto& emdata) {
      const auto Ex = mdspan_t{emdata.Ex.data(), {std::extents{Nx - 1, Ny, Nz}, std::array{Ny * Nz, Nz, 1zu}}};
      const auto Ey = mdspan_t{emdata.Ey.data(), {std::extents{Nx, Ny - 1, Nz}, std::array{(Ny - 1) * Nz, Nz, 1zu}}};
      const auto Ez = mdspan_t{emdata.Ez.data(), {std::extents{Nx, Ny, Nz - 1}, std::array{Ny * (Nz - 1), Nz - 1, 1zu}}};
      const auto Hx = mdspan_t{emdata.Hx.data(), {std::extents{Nx, Ny - 1, Nz - 1}, std::array{(Ny - 1) * (Nz - 1), Nz - 1, 1zu}}};
      const auto Hy = mdspan_t{emdata.Hy.data(), {std::extents{Nx - 1, Ny, Nz - 1}, std::array{Ny * (Nz - 1), Nz - 1, 1zu}}};
      const auto Hz = mdspan_t{emdata.Hz.data(), {std::extents{Nx - 1, Ny - 1, Nz}, std::array{(Ny - 1) * Nz, Nz, 1zu}}};
      
      constexpr auto Chxy = dt / (constants::mu0<double> * dz);
      constexpr auto Chxz = dt / (constants::mu0<double> * dy);
      constexpr auto Chyz = dt / (constants::mu0<double> * dx);
      
      FieldIntegrator<Derivative::FDz, Derivative::FDy>::apply(Hx, Ey, Ez, empty{}, Chxy, Chxz, 0.0);
      FieldIntegrator<Derivative::FDx, Derivative::FDz>::apply(Hy, Ez, Ex, empty{}, Chyz, Chxy, 0.0);
      FieldIntegrator<Derivative::FDy, Derivative::FDx>::apply(Hz, Ex, Ey, empty{}, Chxz, Chyz, 0.0);
   }

   // static void particle_correction(auto& emdata) {
   //    constexpr auto   Bx =   mdspan_t{&emdata.Bx[0, 0, 0], {  hx_shape{}, stride_ijk}};
   //    constexpr auto Eydz = mdspan_t{&emdata.Ey[0, 0, 0], {eydz_shape{}, stride_ijk}};
   //    constexpr auto Ezdy = mdspan_t{&emdata.Ez[0, 0, 0], {ezdy_shape{}, stride_ikj}};
   //
   //    constexpr auto   By =   mdspan_t{&emdata.By[0, 0, 0], {  hy_shape{}, stride_ijk}};
   //    constexpr auto Ezdx = mdspan_t{&emdata.Ez[0, 0, 0], {ezdx_shape{}, stride_kji}};
   //    constexpr auto Exdz = mdspan_t{&emdata.Ex[0, 0, 0], {exdz_shape{}, stride_ijk}};
   //
   //    constexpr auto   Bz =   mdspan_t{&emdata.Bz[0, 0, 0], {  hz_shape{}, stride_ijk}};
   //    constexpr auto Exdy = mdspan_t{&emdata.Ex[0, 0, 0], {exdy_shape{}, stride_ikj}};
   //    constexpr auto Eydx = mdspan_t{&emdata.Ey[0, 0, 0], {eydx_shape{}, stride_kji}};
   //
   //    constexpr auto Chxy2 = 0.5 * dt / (constants::mu0<double> * dz);
   //    constexpr auto Chxz2 = 0.5 * dt / (constants::mu0<double> * dy);
   //    constexpr auto Chyz2 = 0.5 * dt / (constants::mu0<double> * dx);
   //
   //    std::ranges::copy(emdata.Hx, emdata.Bx.begin());
   //    std::ranges::copy(emdata.Hy, emdata.By.begin());
   //    std::ranges::copy(emdata.Hz, emdata.Bz.begin());
   //
   //    Update::apply(Bx, Eydz, Ezdy, empty{}, Chxy2, Chxz2, 0.0);
   //    Update::apply(By, Ezdx, Exdz, empty{}, Chyz2, Chxy2, 0.0);
   //    Update::apply(Bz, Exdy, Eydx, empty{}, Chxz2, Chyz2, 0.0);
   // }

   // static void updateEBCs(auto& emdata) {
   //    constexpr auto Cdx = dt / (constants::eps0<double> * dx);
   //    constexpr auto Cdy = dt / (constants::eps0<double> * dy);
   //    constexpr auto Cdz = dt / (constants::eps0<double> * dz);
   //
   //    X0Ey::apply(Eyx0_psi, Eyx0_b, Eyx0_c, Eyx0, Hzx0_bck, Cdx);
   //    X0Ez::apply(Ezx0_psi, Ezx0_b, Ezx0_c, Ezx0, Hyx0_bck, Cdx);
   //
   //    X1Ey::apply(Eyx1_psi, Eyx1_b, Eyx1_c, Eyx1, Hzx1_bck, Cdx);
   //    X1Ez::apply(Ezx1_psi, Ezx1_b, Ezx1_c, Ezx1, Hyx1_bck, Cdx);
   //
   //    Y0Ex::apply(Exy0_psi, Exy0_b, Exy0_c, Exy0, Hzy0_bck, Cdy);
   //    Y1Ex::apply(Exy1_psi, Exy1_b, Exy1_c, Exy1, Hzy1_bck, Cdy);
   //
   //    Y0Ez::apply(Ezy0_psi, Ezy0_b, Ezy0_c, Ezy0, Hzy0_bck, Cdy);
   //    Y1Ez::apply(Ezy1_psi, Ezy1_b, Ezy1_c, Ezy1, Hzy1_bck, Cdy);
   //
   //    Z0Ex::apply(Exz0_psi, Exz0_b, Exz0_c, Exz0, Hyz0_bck, Cdz);
   //    Z1Ex::apply(Exz1_psi, Exz1_b, Exz1_c, Exz1, Hyz1_bck, Cdz);
   //
   //    Z0Ey::apply(Eyz0_psi, Eyz0_b, Eyz0_c, Eyz0, Hxz0_bck, Cdz);
   //    Z1Ey::apply(Eyz1_psi, Eyz1_b, Eyz1_c, Eyz1, Hxz1_bck, Cdz);
   // }
   //
   // static void updateHBCs(auto& emdata) {
   //    // constexpr auto Chxy = dt / (constants::mu0<double> * dz);
   //    // constexpr auto Chxz = dt / (constants::mu0<double> * dy);
   //    // constexpr auto Chyz = dt / (constants::mu0<double> * dx);
   //    //
   //    // X0Hy::apply(Hyx0, Ezx0, Chyz, emdata.Hyx0);
   //    // X1Hy::apply(Hyx1, Ezx1, Chyz, emdata.Hyx1);
   //    //
   //    // X0Hz::apply(Hzx0, Eyx0, Chyz, emdata.Hzx0);
   //    // X1Hz::apply(Hzx1, Eyx1, Chyz, emdata.Hzx1);
   //    //
   //    // Y0Hx::apply(Hxy0, Ezy0, Chxz, emdata.Hxy0);
   //    // Y1Hx::apply(Hxy1, Ezy1, Chxz, emdata.Hxy1);
   //    //
   //    // Y0Hz::apply(Hzy0, Exy0, Chxz, emdata.Hzy0);
   //    // Y1Hz::apply(Hzy1, Exy1, Chxz, emdata.Hzy1);
   //    //
   //    // Z0Hx::apply(Hxz0, Eyz0, Chxy, emdata.Hxz0);
   //    // Z1Hx::apply(Hxz1, Eyz1, Chxy, emdata.Hxz1);
   //    //
   //    // Z0Hy::apply(Hyz0, Exz0, Chxy, emdata.Hyz0);
   //    // Z1Hy::apply(Hyz1, Exz1, Chxy, emdata.Hyz1);
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
