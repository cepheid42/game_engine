#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "bc_functors.hpp"
// #include "em_data.hpp"
// #include "diff_operators.hpp"

namespace tf::electromagnetics {

template<typename CurlA, typename CurlB>
struct FieldIntegrator {
   using offset_t = std::array<std::size_t, 6>;

   static void operator()(auto& f, const auto& d1, const auto& d2, const auto& src,
                          const auto c_d1, const auto c_d2, const auto c_src)
   {
      #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = 0; i < f.extent(0); ++i) {
         for (std::size_t j = 0; j < f.extent(1); ++j) {
            for (std::size_t k = 0; k < f.extent(2); ++k) {
               const auto diff1   = c_d1 * CurlA::apply(d1, i, j, k);
               const auto diff2   = c_d2 * CurlB::apply(d2, i, j, k);
               const auto current = c_src * src[i, j, k];
               f[i, j, k] = f[i, j, k] + (diff1 - diff2) - current;           
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
}; // end struct FieldIntegrator


struct EMSolver {
   using ExUpdate = FieldIntegrator<backward_dy<Ny>, backward_dz<Nz>>;
   using EyUpdate = FieldIntegrator<backward_dz<Nz>, backward_dx<Nx>>;
   using EzUpdate = FieldIntegrator<backward_dx<Nx>, backward_dy<Ny>>;
   using HxUpdate = FieldIntegrator<forward_dz<Nz>, forward_dy<Ny>>;
   using HyUpdate = FieldIntegrator<forward_dx<Nx>, forward_dz<Nz>>;
   using HzUpdate = FieldIntegrator<forward_dy<Ny>, forward_dx<Nx>>;

   using X0Ey = BCIntegrator<x0bc_t, backward_dx<Nx>, true, true>;
   using X0Ez = BCIntegrator<x0bc_t, backward_dx<Nx>, true, false>;
   using X0Hy = BCIntegrator<x0bc_t, forward_dx<Nx>, true, false>;
   using X0Hz = BCIntegrator<x0bc_t, forward_dx<Nx>, true, true>;
   
   using X1Ey = BCIntegrator<x1bc_t, backward_dx<Nx>, false, true>;
   using X1Ez = BCIntegrator<x1bc_t, backward_dx<Nx>, false, false>;
   using X1Hy = BCIntegrator<x1bc_t, forward_dx<Nx>, false, false>;
   using X1Hz = BCIntegrator<x1bc_t, forward_dx<Nx>, false, true>;
   
   using Y0Ex = BCIntegrator<y0bc_t, backward_dy<Ny>, true, false>;
   using Y0Ez = BCIntegrator<y0bc_t, backward_dy<Ny>, true, true>;
   using Y0Hx = BCIntegrator<y0bc_t, forward_dy<Ny>, true, true>;
   using Y0Hz = BCIntegrator<y0bc_t, forward_dy<Ny>, true, false>;
   
   using Y1Ex = BCIntegrator<y1bc_t, backward_dy<Ny>, false, false>;
   using Y1Ez = BCIntegrator<y1bc_t, backward_dy<Ny>, false, true>;
   using Y1Hx = BCIntegrator<y1bc_t, forward_dy<Ny>, false, true>;
   using Y1Hz = BCIntegrator<y1bc_t, forward_dy<Ny>, false, false>;
   
   using Z0Ex = BCIntegrator<z0bc_t, backward_dz<Nz>, true, true>;
   using Z0Ey = BCIntegrator<z0bc_t, backward_dz<Nz>, true, false>;
   using Z0Hx = BCIntegrator<z0bc_t, forward_dz<Nz>, true, false>;
   using Z0Hy = BCIntegrator<z0bc_t, forward_dz<Nz>, true, true>;
   
   using Z1Ex = BCIntegrator<z1bc_t, backward_dz<Nz>, true, true>;
   using Z1Ey = BCIntegrator<z1bc_t, backward_dz<Nz>, true, false>;
   using Z1Hx = BCIntegrator<z1bc_t, forward_dz<Nz>, true, false>;
   using Z1Hy = BCIntegrator<z1bc_t, forward_dz<Nz>, true, true>;

   struct empty {static auto operator[](auto, auto, auto) { return 0.0; }};
   
   static void advance(auto& emdata, const auto t) requires(em_enabled) {
      updateH(emdata);
      // updateHBCs();
      // updateJBCs();
      // apply_srcs(t);
      updateE(emdata);
      // updateEBCs();
      // particle_correction(); // for the particles and shit
      // zero_currents();       // also for the particles, don't need last week's currents
   }

   static void advance(auto&, const auto) requires (!em_enabled) {}

   static void updateE(auto& emdata) {
      // todo: &emdata.XX[...] are not mdspans, so this will be unhappy
      const auto Ex = mdspan_t{&emdata.Ex[0, 1, 1], std::extents{Nx - 1, Ny - 1, Nz - 1}};
      const auto Ey = mdspan_t{&emdata.Ey[1, 0, 1], std::extents{Nx - 1, Ny - 1, Nz - 1}};
      const auto Ez = mdspan_t{&emdata.Ez[1, 1, 0], std::extents{Nx - 1, Ny - 1, Nz - 1}};
      const auto Jx = mdspan_t{&emdata.Jx[0, 1, 1], std::extents{Nx - 1, Ny - 1, Nz - 1}};
      const auto Jy = mdspan_t{&emdata.Jy[1, 0, 1], std::extents{Nx - 1, Ny - 1, Nz - 1}};
      const auto Jz = mdspan_t{&emdata.Jz[1, 1, 0], std::extents{Nx - 1, Ny - 1, Nz - 1}};
      const auto Hx = mdspan_t{&emdata.Hx[0, 0, 0], std::extents{Nx, Ny - 1, Nz - 1}};
      const auto Hy = mdspan_t{&emdata.Hy[0, 0, 0], std::extents{Nx - 1, Ny, Nz - 1}};
      const auto Hz = mdspan_t{&emdata.Hz[0, 0, 0], std::extents{Nx - 1, Ny - 1, Nz}};
      constexpr auto Cexy = dt / (constants::eps0<double> * dz);
      constexpr auto Cexz = dt / (constants::eps0<double> * dy);
      constexpr auto Ceyz = dt / (constants::eps0<double> * dx);
      constexpr auto Cj = dt / constants::eps0<double>;
      
      ExUpdate(Ex, Hz, Hy, Jx, Ceyz, Cexy, Cj);
      EyUpdate(Ey, Hx, Hz, Jy, Cexy, Ceyz, Cj);
      EzUpdate(Ez, Hy, Hx, Jz, Ceyz, Cexz, Cj);
   }

   static void updateH(auto& emdata) {
      const auto Hx = mdspan_t{&emdata.Hx[0, 0, 0], std::extents{Nx, Ny - 1, Nz - 1}};
      const auto Hy = mdspan_t{&emdata.Hy[0, 0, 0], std::extents{Nx - 1, Ny, Nz - 1}};
      const auto Hz = mdspan_t{&emdata.Hz[0, 0, 0], std::extents{Nx - 1, Ny - 1, Nz}};
      const auto Ex = mdspan_t{&emdata.Ex[0, 0, 0], std::extents{Nx - 1, Ny, Nz}};
      const auto Ey = mdspan_t{&emdata.Ey[0, 0, 0], std::extents{Nx, Ny - 1, Nz}};
      const auto Ez = mdspan_t{&emdata.Ez[0, 0, 0], std::extents{Nx, Ny, Nz - 1}};
      constexpr auto Chxy = dt / (constants::mu0<double> * dz);
      constexpr auto Chxz = dt / (constants::mu0<double> * dy);
      constexpr auto Chyz = dt / (constants::mu0<double> * dx);  
      
      HxUpdate(Hx, Ey, Ez, empty{}, Chxy, Chxz, 0.0);
      HyUpdate(Hy, Ez, Ex, empty{}, Chyz, Chxy, 0.0);
      HzUpdate(Hz, Ex, Ey, empty{}, Chxz, Chyz, 0.0);
   }

   void particle_correction(auto& emdata) {
      const auto Bx = mdspan_t{&emdata.Bx[0, 0, 0], std::extents{Nx, Ny - 1, Nz - 1}};
      const auto By = mdspan_t{&emdata.By[0, 0, 0], std::extents{Nx - 1, Ny, Nz - 1}};
      const auto Bz = mdspan_t{&emdata.Bz[0, 0, 0], std::extents{Nx - 1, Ny - 1, Nz}};
      const auto Ex = mdspan_t{&emdata.Ex[0, 0, 0], std::extents{Nx - 1, Ny, Nz}};
      const auto Ey = mdspan_t{&emdata.Ey[0, 0, 0], std::extents{Nx, Ny - 1, Nz}};
      const auto Ez = mdspan_t{&emdata.Ez[0, 0, 0], std::extents{Nx, Ny, Nz - 1}};
      constexpr auto Chxy2 = 0.5 * dt / (constants::mu0<double> * dz);
      constexpr auto Chxz2 = 0.5 * dt / (constants::mu0<double> * dy);
      constexpr auto Chyz2 = 0.5 * dt / (constants::mu0<double> * dx);
      std::ranges::copy(emdata.Hx, emdata.Bx.begin());
      std::ranges::copy(emdata.Hy, emdata.By.begin());
      std::ranges::copy(emdata.Hz, emdata.Bz.begin());
   
      HxUpdate(Bx, Ey, Ez, empty{}, Chxy2, Chxz2, 0.0);
      HyUpdate(By, Ez, Ex, empty{}, Chyz2, Chxy2, 0.0);
      HzUpdate(Bz, Ex, Ey, empty{}, Chxz2, Chyz2, 0.0);
   }
   
   void updateEBCs(auto& emdata) {
      const auto Eyx0 = mdspan_t{&emdata.Ey[0, 0, 0], std::extents{PMLDepth, Ny - 1, Nz}};
      const auto Eyx1 = mdspan_t{&emdata.Ey[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny - 1, Nz}};

      const auto Hzx0 = mdspan_t{&emdata.Hz[0, 0, 0], std::extents{PMLDepth, Ny - 1, Nz}};
      const auto Hzx1 = mdspan_t{&emdata.Hz[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny - 1, Nz}};

      constexpr auto Cexy = dt / (constants::eps0<double> * dz);
      constexpr auto Cexz = dt / (constants::eps0<double> * dy);
      constexpr auto Ceyz = dt / (constants::eps0<double> * dx);

      X0Ey::apply(Eyx0, Hzx0, Ceyz, emdata.x0_Ey);
      X1Ey::apply(Eyx1, Hzx1, Ceyz, emdata.x1_Ey);
   
      X0Ez::apply(emdata.Ez, emdata.Hy, emdata.Cezhy, emdata.x0.Ez);
      X1Ez::apply(emdata.Ez, emdata.Hy, emdata.Cezhy, emdata.x1.Ez);
   
      Y0Ex::apply(emdata.Ex, emdata.Hz, emdata.Cexhz, emdata.y0.Ex);
      Y1Ex::apply(emdata.Ex, emdata.Hz, emdata.Cexhz, emdata.y1.Ex);
   
      Y0Ez::apply(emdata.Ez, emdata.Hx, emdata.Cezhx, emdata.y0.Ez);
      Y1Ez::apply(emdata.Ez, emdata.Hx, emdata.Cezhx, emdata.y1.Ez);
   
      Z0Ex::apply(emdata.Ex, emdata.Hy, emdata.Cexhy, emdata.z0.Ex);
      Z1Ex::apply(emdata.Ex, emdata.Hy, emdata.Cexhy, emdata.z1.Ex);
   
      Z0Ey::apply(emdata.Ey, emdata.Hx, emdata.Ceyhx, emdata.z0.Ey);
      Z1Ey::apply(emdata.Ey, emdata.Hx, emdata.Ceyhx, emdata.z1.Ey);
   }
   
   void updateHBCs(auto& emdata) {
      X0Hy::apply(emdata.Hy, emdata.Ez, emdata.Chyez, emdata.x0_Hy);
      X1Hy::apply(emdata.Hy, emdata.Ez, emdata.Chyez, emdata.x1_Hy);
   
      X0Hz::apply(emdata.Hz, emdata.Ey, emdata.Chzey, emdata.x0_Hz);
      X1Hz::apply(emdata.Hz, emdata.Ey, emdata.Chzey, emdata.x1_Hz);
   
      Y0Hx::apply(emdata.Hx, emdata.Ez, emdata.Chxez, emdata.y0_Hx);
      Y1Hx::apply(emdata.Hx, emdata.Ez, emdata.Chxez, emdata.y1_Hx);
   
      Y0Hz::apply(emdata.Hz, emdata.Ex, emdata.Chzex, emdata.y0_Hz);
      Y1Hz::apply(emdata.Hz, emdata.Ex, emdata.Chzex, emdata.y1_Hz);
   
      Z0Hx::apply(emdata.Hx, emdata.Ey, emdata.Chxey, emdata.z0_Hx);
      Z1Hx::apply(emdata.Hx, emdata.Ey, emdata.Chxey, emdata.z1_Hx);
   
      Z0Hy::apply(emdata.Hy, emdata.Ex, emdata.Chyex, emdata.z0_Hy);
      Z1Hy::apply(emdata.Hy, emdata.Ex, emdata.Chyex, emdata.z1_Hy);
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
