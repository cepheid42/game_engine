#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

// #include "array_utils.hpp"
#include "constants.hpp"
// #include "bc_functors.hpp"
#include "em_params.hpp"
#include "em_sources.hpp"
#include "program_params.hpp"

// #include <print>
// #include <cassert>

namespace tf::electromagnetics {
template<Derivative D>
static auto diff(const auto& d, const auto i, const auto j, const auto k) {
   if      constexpr (D == Derivative::Dx) { return d[i + 1, j, k] - d[i, j, k]; }
   else if constexpr (D == Derivative::Dy) { return d[i, j + 1, k] - d[i, j, k]; }
   else if constexpr (D == Derivative::Dz) { return d[i, j, k + 1] - d[i, j, k]; }
}

template<Derivative D1, Derivative D2>
struct FieldIntegrator {
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

template<Derivative D>
struct BCIntegrator {
   static constexpr void apply(const auto&, const auto&, const auto&, const auto&) {}

   static void apply(auto& psi, const auto& b, const auto& c, auto& f, const auto& d, const auto& cf)
   {
      for (auto i = 0zu; i < psi.extent(0); ++i) {
         for (auto j = 0zu; j < psi.extent(1); ++j) {
            for (auto k = 0zu; k < psi.extent(2); ++k) {
               psi[i, j, k] = b[k] * psi[i, j, k] + c[k] * diff<D>(d, i, j, k);
               f[i, j, k] += cf * psi[i, j, k]; // cf comes with -1 baked in for += or -=
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
}; // end struct BCIntegrator

struct EMSolver {
   struct empty {
      static constexpr auto operator[](const auto, const auto, const auto) { return 0.0; }
   };
   //                                                    isLo, Negate
   // using X0Ey = BCIntegrator<x0bc_t, backward_dx<Nx>, true, true>;
   // using X0Ez = BCIntegrator<x0bc_t, backward_dx<Nx>, true, false>;
   // using X0Hy = BCIntegrator<x0bc_t, forward_dx<Nx>, true, false>;
   // using X0Hz = BCIntegrator<x0bc_t, forward_dx<Nx>, true, true>;
   //
   // using X1Ey = BCIntegrator<x1bc_t, backward_dx<Nx>, false, true>;
   // using X1Ez = BCIntegrator<x1bc_t, backward_dx<Nx>, false, false>;
   // using X1Hy = BCIntegrator<x1bc_t, forward_dx<Nx>, false, false>;
   // using X1Hz = BCIntegrator<x1bc_t, forward_dx<Nx>, false, true>;
   //
   // using Y0Ex = BCIntegrator<y0bc_t, backward_dy<Ny>, true, false>;
   // using Y0Ez = BCIntegrator<y0bc_t, backward_dy<Ny>, true, true>;
   // using Y0Hx = BCIntegrator<y0bc_t, forward_dy<Ny>, true, true>;
   // using Y0Hz = BCIntegrator<y0bc_t, forward_dy<Ny>, true, false>;
   //
   // using Y1Ex = BCIntegrator<y1bc_t, backward_dy<Ny>, false, false>;
   // using Y1Ez = BCIntegrator<y1bc_t, backward_dy<Ny>, false, true>;
   // using Y1Hx = BCIntegrator<y1bc_t, forward_dy<Ny>, false, true>;
   // using Y1Hz = BCIntegrator<y1bc_t, forward_dy<Ny>, false, false>;
   //
   // using Z0Ex = BCIntegrator<z0bc_t, backward_dz<Nz>, true, true>;
   // using Z0Ey = BCIntegrator<z0bc_t, backward_dz<Nz>, true, false>;
   // using Z0Hx = BCIntegrator<z0bc_t, forward_dz<Nz>, true, false>;
   // using Z0Hy = BCIntegrator<z0bc_t, forward_dz<Nz>, true, true>;
   //
   // using Z1Ex = BCIntegrator<z1bc_t, backward_dz<Nz>, true, true>;
   // using Z1Ey = BCIntegrator<z1bc_t, backward_dz<Nz>, true, false>;
   // using Z1Hx = BCIntegrator<z1bc_t, forward_dz<Nz>, true, false>;
   // using Z1Hy = BCIntegrator<z1bc_t, forward_dz<Nz>, true, true>;

   
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
      const auto Exf = mdspan_t{emdata.Ex.data(), {ex_ext, ex_stride}};
      const auto Eyf = mdspan_t{emdata.Ey.data(), {ey_ext, ey_stride}};
      const auto Ezf = mdspan_t{emdata.Ez.data(), {ez_ext, ez_stride}};
      const auto Jxf = mdspan_t{emdata.Jx.data(), {ex_ext, ex_stride}};
      const auto Jyf = mdspan_t{emdata.Jy.data(), {ey_ext, ey_stride}};
      const auto Jzf = mdspan_t{emdata.Jz.data(), {ez_ext, ez_stride}};
      const auto Hxf = mdspan_t{emdata.Hx.data(), {hx_ext, hx_stride}};
      const auto Hyf = mdspan_t{emdata.Hy.data(), {hy_ext, hy_stride}};
      const auto Hzf = mdspan_t{emdata.Hz.data(), {hz_ext, hz_stride}};

      const auto   Ex = mdspan_t{&Exf[0, 1, 1], {ex_update_ext, ex_stride}};
      const auto   Jx = mdspan_t{&Jxf[0, 1, 1], {ex_update_ext, ex_stride}};
      const auto Hzdy = mdspan_t{&Hzf[0, 0, 1], {ex_update_ext, hz_stride}};
      const auto Hydz = mdspan_t{&Hyf[0, 1, 0], {ex_update_ext, hy_stride}};

      const auto   Ey = mdspan_t{&Eyf[1, 0, 1], {ey_update_ext, ey_stride}};
      const auto   Jy = mdspan_t{&Jyf[1, 0, 1], {ey_update_ext, ey_stride}};
      const auto Hxdz = mdspan_t{&Hxf[1, 0, 0], {ey_update_ext, hx_stride}};
      const auto Hzdx = mdspan_t{&Hzf[0, 0, 1], {ey_update_ext, hz_stride}};

      const auto   Ez = mdspan_t{&Ezf[1, 1, 0], {ez_update_ext, ez_stride}};
      const auto   Jz = mdspan_t{&Jzf[1, 1, 0], {ez_update_ext, ez_stride}};
      const auto Hydx = mdspan_t{&Hyf[0, 1, 0], {ez_update_ext, hy_stride}};
      const auto Hxdy = mdspan_t{&Hxf[1, 0, 0], {ez_update_ext, hx_stride}};

      constexpr auto Cexy = dt / (constants::eps0<double> * dz);
      constexpr auto Cexz = dt / (constants::eps0<double> * dy);
      constexpr auto Ceyz = dt / (constants::eps0<double> * dx);
      constexpr auto Cj = dt / constants::eps0<double>;

      FieldIntegrator<Derivative::Dy, Derivative::Dz>::apply(Ex, Hzdy, Hydz, Jx, Cexz, Ceyz, Cj);
      FieldIntegrator<Derivative::Dz, Derivative::Dx>::apply(Ey, Hxdz, Hzdx, Jy, Cexy, Ceyz, Cj);
      FieldIntegrator<Derivative::Dx, Derivative::Dy>::apply(Ez, Hydx, Hxdy, Jz, Ceyz, Cexz, Cj);
   }

   static void updateH(auto& emdata) {
      const auto Ex = mdspan_t{emdata.Ex.data(), {ex_ext, ex_stride}};
      const auto Ey = mdspan_t{emdata.Ey.data(), {ey_ext, ey_stride}};
      const auto Ez = mdspan_t{emdata.Ez.data(), {ez_ext, ez_stride}};
      const auto Hx = mdspan_t{emdata.Hx.data(), {hx_ext, hx_stride}};
      const auto Hy = mdspan_t{emdata.Hy.data(), {hy_ext, hy_stride}};
      const auto Hz = mdspan_t{emdata.Hz.data(), {hz_ext, hz_stride}};

      constexpr auto Chxy = dt / (constants::mu0<double> * dz);
      constexpr auto Chxz = dt / (constants::mu0<double> * dy);
      constexpr auto Chyz = dt / (constants::mu0<double> * dx);
      
      FieldIntegrator<Derivative::Dz, Derivative::Dy>::apply(Hx, Ey, Ez, empty{}, Chxy, Chxz, 0.0);
      FieldIntegrator<Derivative::Dx, Derivative::Dz>::apply(Hy, Ez, Ex, empty{}, Chyz, Chxy, 0.0);
      FieldIntegrator<Derivative::Dy, Derivative::Dx>::apply(Hz, Ex, Ey, empty{}, Chxz, Chyz, 0.0);
   }

   static void particle_correction(auto& emdata) {
      const auto Ex = mdspan_t{emdata.Ex.data(), {ex_ext, ex_stride}};
      const auto Ey = mdspan_t{emdata.Ey.data(), {ey_ext, ey_stride}};
      const auto Ez = mdspan_t{emdata.Ez.data(), {ez_ext, ez_stride}};
      const auto Bx = mdspan_t{emdata.Hx.data(), {hx_ext, hx_stride}};
      const auto By = mdspan_t{emdata.Hy.data(), {hy_ext, hy_stride}};
      const auto Bz = mdspan_t{emdata.Hz.data(), {hz_ext, hz_stride}};

      constexpr auto Chxy2 = 0.5 * dt / (constants::mu0<double> * dz);
      constexpr auto Chxz2 = 0.5 * dt / (constants::mu0<double> * dy);
      constexpr auto Chyz2 = 0.5 * dt / (constants::mu0<double> * dx);

      std::ranges::copy(emdata.Hx, emdata.Bx.begin());
      std::ranges::copy(emdata.Hy, emdata.By.begin());
      std::ranges::copy(emdata.Hz, emdata.Bz.begin());

      FieldIntegrator<Derivative::Dz, Derivative::Dy>::apply(Bx, Ey, Ez, empty{}, Chxy2, Chxz2, 0.0);
      FieldIntegrator<Derivative::Dx, Derivative::Dz>::apply(By, Ez, Ex, empty{}, Chyz2, Chxy2, 0.0);
      FieldIntegrator<Derivative::Dy, Derivative::Dx>::apply(Bz, Ex, Ey, empty{}, Chxz2, Chyz2, 0.0);
   }

   static void updateEBCs(auto& emdata) {
      constexpr auto Cdx = dt / (constants::eps0<double> * dx);
      // constexpr auto Cdy = dt / (constants::eps0<double> * dy);
      // constexpr auto Cdz = dt / (constants::eps0<double> * dz);

      // const auto Exf = mdspan_t{emdata.Ex.data(), {ex_ext, ex_stride}};
      const auto Eyf = mdspan_t{emdata.Ey.data(), {ey_ext, ey_stride}};
      const auto Ezf = mdspan_t{emdata.Ez.data(), {ez_ext, ez_stride}};
      // const auto Hxf = mdspan_t{emdata.Hx.data(), {hx_ext, hx_stride}};
      const auto Hyf = mdspan_t{emdata.Hy.data(), {hy_ext, hy_stride}};
      const auto Hzf = mdspan_t{emdata.Hz.data(), {hz_ext, hz_stride}};

      constexpr auto eyx0_ext = std::extents{BCDepth, Ny - 1, Nz};
      constexpr auto hzx0_ext = std::extents{BCDepth, Ny - 1, Nz};
      const auto eyx0_psi = mdspan_t{emdata.Eyx0.psi.data(), {eyx0_ext, ey_stride}};
      const auto eyx0_b = std::span{emdata.Eyx0.b.data()};
      const auto eyx0_c = std::span{emdata.Eyx0.c.data()};
      const auto eyx0 = mdspan_t{&Eyf[0, 0, 0], {eyx0_ext, ey_stride}};
      const auto hzx0 = mdspan_t{&Hzf[0, 0, 0], {hzx0_ext, hz_stride}};
      const auto eyx1_psi = mdspan_t{emdata.Eyx1.data(), {eyx0_ext, ey_stride}};
      const auto eyx1_b = std::span{emdata.Eyx1.b.data()};
      const auto eyx1_c = std::span{emdata.Eyx1.c.data()};
      const auto eyx1 = mdspan_t{&Eyf[Nx - BCDepth, 0, 0], {eyx0_ext, ey_stride}};
      const auto hzx1 = mdspan_t{&Hzf[Nx - BCDepth - 1, 0, 0], {hzx0_ext, hz_stride}}; // todo: need the -1?
      BCIntegrator<Derivative::Dx>::apply(eyx0_psi, eyx0_b, eyx0_c, eyx0, hzx0, -Cdx);
      BCIntegrator<Derivative::Dx>::apply(eyx1_psi, eyx1_b, eyx1_c, eyx1, hzx1, -Cdx);


      constexpr auto ezx0_ext = std::extents{BCDepth, Ny, Nz - 1};
      constexpr auto hyx0_ext = std::extents{BCDepth, Ny, Nz - 1};
      const auto ezx0_psi = mdspan_t{emdata.Eyx0.psi.data(), {ezx0_ext, ez_stride}};
      const auto ezx0_b = std::span{emdata.Eyx0.b.data()};
      const auto ezx0_c = std::span{emdata.Eyx0.c.data()};
      const auto ezx0 = mdspan_t{&Ezf[0, 0, 0], {ezx0_ext, ez_stride}};
      const auto hyx0 = mdspan_t{&Hyf[0, 0, 0], {hyx0_ext, hy_stride}};
      const auto ezx1_psi = mdspan_t{emdata.Eyx1.psi.data(), {ezx0_ext, ez_stride}};
      const auto ezx1_b = std::span{emdata.Eyx1.b.data()};
      const auto ezx1_c = std::span{emdata.Eyx1.c.data()};
      const auto ezx1 = mdspan_t{&Ezf[Nx - BCDepth, 0, 0], {ezx0_ext, ez_stride}};
      const auto hyx1 = mdspan_t{&Hyf[Nx - BCDepth - 1, 0, 0], {hyx0_ext, hy_stride}}; // todo: need the -1?
      BCIntegrator<Derivative::Dx>::apply(ezx0_psi, ezx0_b, ezx0_c, ezx0, hyx0, Cdx);
      BCIntegrator<Derivative::Dx>::apply(ezx1_psi, ezx1_b, ezx1_c, ezx1, hyx1, Cdx);

      // BCIntegrator<Derivative::Dy>::apply(Exy0_psi, Exy0_b, Exy0_c, Ex, Hzy0_bck, Cdy);
      // BCIntegrator<Derivative::Dy>::apply(Exy1_psi, Exy1_b, Exy1_c, Ex, Hzy1_bck, Cdy);
      //
      // BCIntegrator<Derivative::Dy>::apply(Ezy0_psi, Ezy0_b, Ezy0_c, Ez, Hzy0_bck, -Cdy);
      // BCIntegrator<Derivative::Dy>::apply(Ezy1_psi, Ezy1_b, Ezy1_c, Ez, Hzy1_bck, -Cdy);
      //
      // BCIntegrator<Derivative::Dz>::apply(Exz0_psi, Exz0_b, Exz0_c, Ex, Hyz0_bck, -Cdz);
      // BCIntegrator<Derivative::Dz>::apply(Exz1_psi, Exz1_b, Exz1_c, Ex, Hyz1_bck, -Cdz);
      //
      // BCIntegrator<Derivative::Dz>::apply(Eyz0_psi, Eyz0_b, Eyz0_c, Ey, Hxz0_bck, Cdz);
      // BCIntegrator<Derivative::Dz>::apply(Eyz1_psi, Eyz1_b, Eyz1_c, Ey, Hxz1_bck, Cdz);
   }

   static void updateHBCs(auto& emdata) {
      constexpr auto Cdx = dt / (constants::mu0<double> * dx);
      // constexpr auto Cdy = dt / (constants::mu0<double> * dy);
      // constexpr auto Cdz = dt / (constants::mu0<double> * dz);

      // const auto Exf = mdspan_t{emdata.Ex.data(), {ex_ext, ex_stride}};
      const auto Eyf = mdspan_t{emdata.Ey.data(), {ey_ext, ey_stride}};
      const auto Ezf = mdspan_t{emdata.Ez.data(), {ez_ext, ez_stride}};
      // const auto Hxf = mdspan_t{emdata.Hx.data(), {hx_ext, hx_stride}};
      const auto Hyf = mdspan_t{emdata.Hy.data(), {hy_ext, hy_stride}};
      const auto Hzf = mdspan_t{emdata.Hz.data(), {hz_ext, hz_stride}};

      constexpr auto hyx0_ext = std::extents{BCDepth, Ny, Nz - 1};
      constexpr auto ezx0_ext = std::extents{BCDepth, Ny, Nz - 1};
      const auto hyx0_psi = mdspan_t{emdata.Hyx0.psi.data(), {hyx0_ext, hy_stride}};
      const auto hyx0_b = std::span{emdata.Hyx0.b.data()};
      const auto hyx0_c = std::span{emdata.Hyx0.c.data()};
      const auto hyx0 = mdspan_t{&Hyf[0, 0, 0], {hyx0_ext, hy_stride}};
      const auto ezx0 = mdspan_t{&Ezf[0, 0, 0], {ezx0_ext, ez_stride}};
      const auto hyx1_psi = mdspan_t{emdata.Hyx1.psi.data(), {hyx0_ext, hy_stride}};
      const auto hyx1_b = std::span{emdata.Hyx1.b.data()};
      const auto hyx1_c = std::span{emdata.Hyx1.c.data()};
      const auto hyx1 = mdspan_t{&Hyf[Nx - BCDepth - 1, 0, 0], {hyx0_ext, hy_stride}}; // todo: need the -1?
      const auto ezx1 = mdspan_t{&Ezf[Nx - BCDepth, 0, 0], {ezx0_ext, ez_stride}};
      BCIntegrator<Derivative::Dx>::apply(hyx0_psi, hyx0_b, hyx0_c, hyx0, ezx0, Cdx);
      BCIntegrator<Derivative::Dx>::apply(hyx1_psi, hyx1_b, hyx1_c, hyx1, ezx1, Cdx);

      constexpr auto hzx0_ext = std::extents{BCDepth, Ny - 1, Nz};
      constexpr auto eyx0_ext = std::extents{BCDepth, Ny - 1, Nz};
      const auto hzx0_psi = mdspan_t{emdata.Hx0.psi.data(), {hzx0_ext, hz_stride}};
      const auto hzx0_b = std::span{emdata.Hzx0.b.data()};
      const auto hzx0_c = std::span{emdata.Hzx0.c.data()};
      const auto hzx0 = mdspan_t{&Hzf[0, 0, 0], {hzx0_ext, hz_stride}};
      const auto eyx0 = mdspan_t{&Eyf[0, 0, 0], {eyx0_ext, ey_stride}};
      const auto hzx1_psi = mdspan_t{emdata.Hzx1.psi.data(), {hzx0_ext, hz_stride}};
      const auto hzx1_b = std::span{emdata.Hzx1.b.data()};
      const auto hzx1_c = std::span{emdata.Hzx1.c.data()};
      const auto hzx1 = mdspan_t{&Hzf[Nx - BCDepth - 1, 0, 0], {hzx0_ext, hz_stride}}; // todo: need the -1?
      const auto eyx1 = mdspan_t{&Eyf[Nx - BCDepth, 0, 0], {eyx0_ext, ey_stride}};
      BCIntegrator<Derivative::Dx>::apply(hzx0_psi, hzx0_b, hzx0_c, hzx0, eyx0, -Cdx);
      BCIntegrator<Derivative::Dx>::apply(hzx1_psi, hzx1_b, hzx1_c, hzx1, eyx1, -Cdx);

      // Y0Hx::apply(Hxy0, Ezy0, Chxz, emdata.Hxy0);
      // Y1Hx::apply(Hxy1, Ezy1, Chxz, emdata.Hxy1);
      //
      // Y0Hz::apply(Hzy0, Exy0, Chxz, emdata.Hzy0);
      // Y1Hz::apply(Hzy1, Exy1, Chxz, emdata.Hzy1);
      //
      // Z0Hx::apply(Hxz0, Eyz0, Chxy, emdata.Hxz0);
      // Z1Hx::apply(Hxz1, Eyz1, Chxy, emdata.Hxz1);
      //
      // Z0Hy::apply(Hyz0, Exz0, Chxy, emdata.Hyz0);
      // Z1Hy::apply(Hyz1, Exz1, Chxy, emdata.Hyz1);
   }
   
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
