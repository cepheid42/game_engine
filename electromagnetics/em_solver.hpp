#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "constants.hpp"
// #include "bc_functors.hpp"
#include "em_params.hpp"
#include "em_sources.hpp"

#include <print>

namespace tf::electromagnetics {

template<Derivative D1, Derivative D2>
struct FieldIntegrator {
   template<Derivative D>
   static auto diff(const auto& d, const auto i, const auto j, const auto k) {
      if      constexpr (D == Derivative::FDx) { return d[i + 1, j, k] - d[i, j, k]; }
      else if constexpr (D == Derivative::FDy) { return d[i, j + 1, k] - d[i, j, k]; }
      else if constexpr (D == Derivative::FDz) { return d[i, j, k + 1] - d[i, j, k]; }
      else if constexpr (D == Derivative::BDx) { return d[i, j, k] - d[i - 1, j, k]; }
      else if constexpr (D == Derivative::BDy) { return d[i, j, k] - d[i, j - 1, k]; }
      else if constexpr (D == Derivative::BDz) { return d[i, j, k] - d[i, j, k - 1]; }
   }

   static void apply(auto& f, const auto& d1, const auto& d2, const auto& src,
                     const auto c_d1, const auto c_d2, const auto c_src, const auto& offsets)
   {
      const auto& [x0, x1, y0, y1, z0, z1] = offsets;
      // #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = x0; i < f.extent(0) - x1; ++i) {
         for (std::size_t j = y0; j < f.extent(1) - y1; ++j) {
            for (std::size_t k = z0; k < f.extent(2) - z1; ++k) {
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
      const auto Ez = mdspan_t{emdata.Ez.data(), Nx, Ny, Nz - 1};

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
      const auto Ex = mdspan_t{emdata.Ex.data(), Nx - 1, Ny, Nz};
      const auto Ey = mdspan_t{emdata.Ey.data(), Nx, Ny - 1, Nz};
      const auto Ez = mdspan_t{emdata.Ez.data(), Nx, Ny, Nz - 1};
      const auto Jx = mdspan_t{emdata.Jx.data(), Nx - 1, Ny, Nz};
      const auto Jy = mdspan_t{emdata.Jy.data(), Nx, Ny - 1, Nz};
      const auto Jz = mdspan_t{emdata.Jz.data(), Nx, Ny, Nz - 1};

      const auto Hx = mdspan_t{emdata.Hx.data(), Nx, Ny - 1, Nz - 1};
      const auto Hy = mdspan_t{emdata.Hy.data(), Nx - 1, Ny, Nz - 1};
      const auto Hz = mdspan_t{emdata.Hz.data(), Nx - 1, Ny - 1, Nz};

      constexpr auto Cexy = dt / (constants::eps0<double> * dz);
      constexpr auto Cexz = dt / (constants::eps0<double> * dy);
      constexpr auto Ceyz = dt / (constants::eps0<double> * dx);
      constexpr auto Cj = dt / constants::eps0<double>;

      FieldIntegrator<Derivative::BDy, Derivative::BDz>::apply(Ex, Hz, Hy, Jx, Cexz, Ceyz, Cj, std::array{0, 0, 1, 1, 1, 1});
      FieldIntegrator<Derivative::BDz, Derivative::BDx>::apply(Ey, Hx, Hz, Jy, Cexy, Ceyz, Cj, std::array{1, 1, 0, 0, 1, 1});
      FieldIntegrator<Derivative::BDx, Derivative::BDy>::apply(Ez, Hy, Hx, Jz, Ceyz, Cexz, Cj, std::array{1, 1, 1, 1, 0, 0});
   }

   static void updateH(auto& emdata) {
      const auto Ex = mdspan_t{emdata.Ex.data(), Nx - 1, Ny, Nz};
      const auto Ey = mdspan_t{emdata.Ey.data(), Nx, Ny - 1, Nz};
      const auto Ez = mdspan_t{emdata.Ez.data(), Nx, Ny, Nz - 1};
      const auto Hx = mdspan_t{emdata.Hx.data(), Nx, Ny - 1, Nz - 1};
      const auto Hy = mdspan_t{emdata.Hy.data(), Nx - 1, Ny, Nz - 1};
      const auto Hz = mdspan_t{emdata.Hz.data(), Nx - 1, Ny - 1, Nz};

      constexpr auto Chxy = dt / (constants::mu0<double> * dz);
      constexpr auto Chxz = dt / (constants::mu0<double> * dy);
      constexpr auto Chyz = dt / (constants::mu0<double> * dx);

      FieldIntegrator<Derivative::FDz, Derivative::FDy>::apply(Hx, Ey, Ez, empty{}, Chxy, Chxz, 0.0, std::array{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
      FieldIntegrator<Derivative::FDx, Derivative::FDz>::apply(Hy, Ez, Ex, empty{}, Chyz, Chxy, 0.0, std::array{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
      FieldIntegrator<Derivative::FDy, Derivative::FDx>::apply(Hz, Ex, Ey, empty{}, Chxz, Chyz, 0.0, std::array{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
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
   //    // X0 EyHz
   //    constexpr auto Eyx0_b = std::span{&emdata.Eyx0.b[0], BCDepth - 1};
   //    constexpr auto Eyx0_c = std::span{&emdata.Eyx0.c[0], BCDepth - 1};
   //    constexpr auto     Eyx0 = mdspan_t{      &emdata.Ey[1, 0, 0], {ey_shape_pml{}, ey_stride_pml}};
   //    constexpr auto Hzx0_bck = mdspan_t{      &emdata.Hz[1, 0, 0], {hzdx_shape_pml{}, hzdx_stride_pml}};
   //    constexpr auto Eyx0_psi = mdspan_t{&emdata.Eyx0.psi[1, 0, 0], {ey_shape_pml{}, ey_stride_pml}};
   //
   //    // X0 EzHy
   //    constexpr auto Ezx0_b = std::span{&emdata.Ezx0.b[0], BCDepth - 1};
   //    constexpr auto Ezx0_c = std::span{&emdata.Ezx0.c[0], BCDepth - 1};
   //    constexpr auto     Ezx0 = mdspan_t{      &emdata.Ez[1, 0, 0], {ez_shape_pml{}  , ez_stride_pml}};
   //    constexpr auto Hyx0_bck = mdspan_t{      &emdata.Hy[1, 0, 0], {hydx_shape_pml{}, hydx_stride_pml}};
   //    constexpr auto Ezx0_psi = mdspan_t{&emdata.Ezx0.psi[1, 0, 0], {ez_shape_pml{}  , ez_stride_pml}};
   //
   //    // X1 EyHz
   //    constexpr auto Eyx1_b = std::span{&emdata.Eyx1.b[0], BCDepth - 1};
   //    constexpr auto Eyx1_c = std::span{&emdata.Eyx1.c[0], BCDepth - 1};
   //    constexpr auto     Eyx1 = mdspan_t{      &emdata.Ey[Nx - BCDepth, 0, 0], {ey_shape_pml{}  , ey_stride_pml}};
   //    constexpr auto Hzx1_bck = mdspan_t{      &emdata.Hz[Nx - BCDepth, 0, 0], {hzdx_shape_pml{}, hzdx_stride_pml}};
   //    constexpr auto Eyx1_psi = mdspan_t{&emdata.Eyx1.psi[Nx - BCDepth, 0, 0], {ey_shape_pml{}  , ey_stride_pml}};
   //
   //    // X1 EzHy
   //    constexpr auto Ezx1_b = std::span{&emdata.Ezx1.b[0], BCDepth - 1};
   //    constexpr auto Ezx1_c = std::span{&emdata.Ezx1.c[0], BCDepth - 1};
   //    constexpr auto     Ezx1 = mdspan_t{      &emdata.Ez[Nx - BCDepth, 0, 0], {ez_shape_pml{}  , ez_stride_pml}};
   //    constexpr auto Hyx1_bck = mdspan_t{      &emdata.Hy[Nx - BCDepth, 0, 0], {hydx_shape_pml{}, hydx_stride_pml}};
   //    constexpr auto Ezx1_psi = mdspan_t{&emdata.Ezx1.psi[Nx - BCDepth, 0, 0], {ez_shape_pml{}  , ez_stride_pml}};
   //
   //    // Y0 ExHz
   //    constexpr auto Exy0_b = std::span{&emdata.Exy0.b[0], BCDepth - 1};
   //    constexpr auto Exy0_c = std::span{&emdata.Exy0.c[0], BCDepth - 1};
   //    constexpr auto     Exy0 = mdspan_t{      &emdata.Ex[0, 1, 0], {ex_shape_pml{}  , ex_stride_pml}};
   //    constexpr auto Hzy0_bck = mdspan_t{      &emdata.Hz[0, 1, 0], {hzdx_shape_pml{}, hzdx_stride_pml}};
   //    constexpr auto Exy0_psi = mdspan_t{&emdata.Exy0.psi[0, 1, 0], {ex_shape_pml{}  , ex_stride_pml}};
   //
   //    // Y0 EzHx
   //    constexpr auto Ezy0_b = std::span{&emdata.Ezy0.b[0], BCDepth - 1};
   //    constexpr auto Ezy0_c = std::span{&emdata.Ezy0.c[0], BCDepth - 1};
   //    constexpr auto     Ezy0 = mdspan_t{      &emdata.Ez[0, 1, 0], {ez_shape_pml{}  , ez_stride_pml}};
   //    constexpr auto Hxy0_bck = mdspan_t{      &emdata.Hx[0, 1, 0], {hxdy_shape_pml{}, hxdy_stride_pml}};
   //    constexpr auto Ezy0_psi = mdspan_t{&emdata.Ezy0.psi[0, 1, 0], {ez_shape_pml{}  , ez_stride_pml}};
   //
   //    // Y1 ExHz
   //    constexpr auto Exy1_b = std::span{&emdata.Exy1.b[0], BCDepth - 1};
   //    constexpr auto Exy1_c = std::span{&emdata.Exy1.c[0], BCDepth - 1};
   //    constexpr auto     Exy1 = mdspan_t{      &emdata.Ex[0, Ny - BCDepth, 0], {ex_shape_pml{}  , ex_stride_pml}};
   //    constexpr auto Hzy1_bck = mdspan_t{      &emdata.Hz[0, Ny - BCDepth, 0], {hzdx_shape_pml{}, hzdx_stride_pml}};
   //    constexpr auto Exy1_psi = mdspan_t{&emdata.Exy1.psi[0, Ny - BCDepth, 0], {ex_shape_pml{}  , ex_stride_pml}};
   //
   //    // Y1 EzHx
   //    constexpr auto Ezy1_b = std::span{&emdata.Ezy1.b[0], BCDepth - 1};
   //    constexpr auto Ezy1_c = std::span{&emdata.Ezy1.c[0], BCDepth - 1};
   //    constexpr auto     Ezy1 = mdspan_t{      &emdata.Ez[0, Ny - BCDepth, 0], {ez_shape_pml{}  , ez_stride_pml}};
   //    constexpr auto Hxy1_bck = mdspan_t{      &emdata.Hx[0, Ny - BCDepth, 0], {hxdy_shape_pml{}, hxdy_stride_pml}};
   //    constexpr auto Ezy1_psi = mdspan_t{&emdata.Ezy1.psi[0, Ny - BCDepth, 0], {ez_shape_pml{}  , ez_stride_pml}};
   //
   //    // Z0 ExHy
   //    constexpr auto Exz0_b = std::span{&emdata.Exz0.b[0], BCDepth - 1};
   //    constexpr auto Exz0_c = std::span{&emdata.Exz0.c[0], BCDepth - 1};
   //    constexpr auto     Exz0 = mdspan_t{      &emdata.Ex[0, 0, 1], {ex_shape_pml{}  , ex_stride_pml}};
   //    constexpr auto Hyz0_bck = mdspan_t{      &emdata.Hy[0, 0, 1], {hydx_shape_pml{}, hydx_stride_pml}};
   //    constexpr auto Exz0_psi = mdspan_t{&emdata.Exz0.psi[0, 0, 1], {ex_shape_pml{}  , ex_stride_pml}};
   //
   //    // z0 EyHx
   //    constexpr auto Eyz0_b = std::span{&emdata.Eyz0.b[0], BCDepth - 1};
   //    constexpr auto Eyz0_c = std::span{&emdata.Eyz0.c[0], BCDepth - 1};
   //    constexpr auto     Eyz0 = mdspan_t{      &emdata.Ey[0, 0, 1], {ey_shape_pml{}  , ey_stride_pml}};
   //    constexpr auto Hxz0_bck = mdspan_t{      &emdata.Hx[0, 0, 1], {hxdy_shape_pml{}, hxdy_stride_pml}};
   //    constexpr auto Eyz0_psi = mdspan_t{&emdata.Eyz0.psi[0, 0, 1], {ey_shape_pml{}  , ey_stride_pml}};
   //
   //    // z1 ExHy
   //    constexpr auto Exz1_b = std::span{&emdata.Exz1.b[0], BCDepth - 1};
   //    constexpr auto Exz1_c = std::span{&emdata.Exz1.c[0], BCDepth - 1};
   //    constexpr auto     Exz1 = mdspan_t{      &emdata.Ex[0, 0, Nz - BCDepth], {ex_shape_pml{}  , ex_stride_pml}};
   //    constexpr auto Hyz1_bck = mdspan_t{      &emdata.Hy[0, 0, Nz - BCDepth], {hzdx_shape_pml{}, hydx_stride_pml}};
   //    constexpr auto Exz1_psi = mdspan_t{&emdata.Exz1.psi[0, 0, Nz - BCDepth], {ex_shape_pml{}  , ex_stride_pml}};
   //
   //    // z1 EyHx
   //    constexpr auto Eyz1_b = std::span{&emdata.Eyz1.b[0], BCDepth - 1};
   //    constexpr auto Eyz1_c = std::span{&emdata.Eyz1.c[0], BCDepth - 1};
   //    constexpr auto     Eyz1 = mdspan_t{      &emdata.Ey[0, 0, Nz - BCDepth], {ey_shape_pml{}  , ey_stride_pml}};
   //    constexpr auto Hxz1_bck = mdspan_t{      &emdata.Hx[0, 0, Nz - BCDepth], {hxdy_shape_pml{}, hxdy_stride_pml}};
   //    constexpr auto Eyz1_psi = mdspan_t{&emdata.Eyz1.psi[0, 0, Nz - BCDepth], {ey_shape_pml{}  , ey_stride_pml}};
   //
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
