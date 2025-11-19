#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "bc_functors.hpp"
// #include "em_data.hpp"
// #include "diff_operators.hpp"

namespace tf::electromagnetics {

template<typename CurlA, typename CurlB>
struct FieldIntegrator {
   using offset_t = std::array<std::size_t, 6>;

   static void apply(auto& f, const auto& d1, const auto& d2, const auto& src,
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
      updateHBCs(emdata);
      // updateJBCs();
      // apply_srcs(t);
      updateE(emdata);
      updateEBCs(emdata);
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
      
      ExUpdate::apply(Ex, Hz, Hy, Jx, Ceyz, Cexy, Cj);
      EyUpdate::apply(Ey, Hx, Hz, Jy, Cexy, Ceyz, Cj);
      EzUpdate::apply(Ez, Hy, Hx, Jz, Ceyz, Cexz, Cj);
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
      
      HxUpdate::apply(Hx, Ey, Ez, empty{}, Chxy, Chxz, 0.0);
      HyUpdate::apply(Hy, Ez, Ex, empty{}, Chyz, Chxy, 0.0);
      HzUpdate::apply(Hz, Ex, Ey, empty{}, Chxz, Chyz, 0.0);
   }

   static void particle_correction(auto& emdata) {
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
   
      HxUpdate::apply(Bx, Ey, Ez, empty{}, Chxy2, Chxz2, 0.0);
      HyUpdate::apply(By, Ez, Ex, empty{}, Chyz2, Chxy2, 0.0);
      HzUpdate::apply(Bz, Ex, Ey, empty{}, Chxz2, Chyz2, 0.0);
   }
   
   static void updateEBCs(auto& emdata) {
      const auto Eyx0 = mdspan_t{&emdata.Ey[0, 0, 0], std::extents{PMLDepth, Ny - 1  , Nz      }};
      const auto Ezx0 = mdspan_t{&emdata.Ez[0, 0, 0], std::extents{PMLDepth, Ny      , Nz - 1  }};
      const auto Exy0 = mdspan_t{&emdata.Ex[0, 0, 0], std::extents{Nx - 1  , PMLDepth, Nz      }};
      const auto Ezy0 = mdspan_t{&emdata.Ez[0, 0, 0], std::extents{Nx      , PMLDepth, Nz - 1  }};
      const auto Exz0 = mdspan_t{&emdata.Ex[0, 0, 0], std::extents{Nx - 1  , Ny      , PMLDepth}};
      const auto Eyz0 = mdspan_t{&emdata.Ey[0, 0, 0], std::extents{Nx      , Ny - 1  , PMLDepth}};

      const auto Eyx1 = mdspan_t{&emdata.Ey[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny - 1  , Nz      }};
      const auto Ezx1 = mdspan_t{&emdata.Ez[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny      , Nz - 1  }};
      const auto Exy1 = mdspan_t{&emdata.Ex[0, Ny - PMLDepth, 0], std::extents{Nx - 1  , PMLDepth, Nz      }};
      const auto Ezy1 = mdspan_t{&emdata.Ez[0, Ny - PMLDepth, 0], std::extents{Nx      , PMLDepth, Nz - 1  }};
      const auto Exz1 = mdspan_t{&emdata.Ex[0, 0, Nz - PMLDepth], std::extents{Nx - 1  , Ny      , PMLDepth}};
      const auto Eyz1 = mdspan_t{&emdata.Ey[0, 0, Nz - PMLDepth], std::extents{Nx      , Ny - 1  , PMLDepth}};

      const auto Hyx0 = mdspan_t{&emdata.Hy[0, 0, 0], std::extents{PMLDepth, Ny      , Nz - 1  }};
      const auto Hzx0 = mdspan_t{&emdata.Hz[0, 0, 0], std::extents{PMLDepth, Ny - 1  , Nz      }};
      const auto Hxy0 = mdspan_t{&emdata.Hx[0, 0, 0], std::extents{Nx      , PMLDepth, Nz - 1  }};
      const auto Hzy0 = mdspan_t{&emdata.Hz[0, 0, 0], std::extents{Nx - 1  , PMLDepth, Nz      }};
      const auto Hxz0 = mdspan_t{&emdata.Hx[0, 0, 0], std::extents{Nx      , Ny - 1  , PMLDepth}};
      const auto Hyz0 = mdspan_t{&emdata.Hy[0, 0, 0], std::extents{Nx - 1  , Ny      , PMLDepth}};
      
      const auto Hyx1 = mdspan_t{&emdata.Hy[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny      , Nz - 1  }};
      const auto Hzx1 = mdspan_t{&emdata.Hz[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny - 1  , Nz      }};
      const auto Hxy1 = mdspan_t{&emdata.Hx[0, Ny - PMLDepth, 0], std::extents{Nx      , PMLDepth, Nz - 1  }};
      const auto Hzy1 = mdspan_t{&emdata.Hz[0, Ny - PMLDepth, 0], std::extents{Nx - 1  , PMLDepth, Nz      }};
      const auto Hxz1 = mdspan_t{&emdata.Hx[0, 0, Nz - PMLDepth], std::extents{Nx      , Ny - 1  , PMLDepth}};
      const auto Hyz1 = mdspan_t{&emdata.Hy[0, 0, Nz - PMLDepth], std::extents{Nx - 1  , Ny      , PMLDepth}};

      constexpr auto Cexy = dt / (constants::eps0<double> * dz);
      constexpr auto Cexz = dt / (constants::eps0<double> * dy);
      constexpr auto Ceyz = dt / (constants::eps0<double> * dx);

      X0Ey::apply(Eyx0, Hzx0, Ceyz, emdata.x0_Ey);
      X1Ey::apply(Eyx1, Hzx1, Ceyz, emdata.x1_Ey);
   
      X0Ez::apply(Ezx0, Hyx0, Ceyz, emdata.x0_Ez);
      X1Ez::apply(Ezx1, Hyx1, Ceyz, emdata.x1_Ez);
   
      Y0Ex::apply(Exy0, Hzy0, Cexz, emdata.y0_Ex);
      Y1Ex::apply(Exy1, Hzy1, Cexz, emdata.y1_Ex);
   
      Y0Ez::apply(Ezy0, Hxy0, Cexz, emdata.y0_Ez);
      Y1Ez::apply(Ezy1, Hxy1, Cexz, emdata.y1_Ez);
   
      Z0Ex::apply(Exz0, Hyz0, Cexy, emdata.z0_Ex);
      Z1Ex::apply(Exz1, Hyz1, Cexy, emdata.z1_Ex);
   
      Z0Ey::apply(Eyz0, Hxz0, Cexy, emdata.z0_Ey);
      Z1Ey::apply(Eyz1, Hxz1, Cexy, emdata.z1_Ey);
   }
   
   static void updateHBCs(auto& emdata) {
      const auto Hyx0 = mdspan_t{&emdata.Hy[0, 0, 0], std::extents{PMLDepth, Ny      , Nz - 1  }};
      const auto Hzx0 = mdspan_t{&emdata.Hz[0, 0, 0], std::extents{PMLDepth, Ny - 1  , Nz      }};
      const auto Hxy0 = mdspan_t{&emdata.Hx[0, 0, 0], std::extents{Nx      , PMLDepth, Nz - 1  }};
      const auto Hzy0 = mdspan_t{&emdata.Hz[0, 0, 0], std::extents{Nx - 1  , PMLDepth, Nz      }};
      const auto Hxz0 = mdspan_t{&emdata.Hx[0, 0, 0], std::extents{Nx      , Ny - 1  , PMLDepth}};
      const auto Hyz0 = mdspan_t{&emdata.Hy[0, 0, 0], std::extents{Nx - 1  , Ny      , PMLDepth}};

      const auto Hyx1 = mdspan_t{&emdata.Hy[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny      , Nz - 1  }};
      const auto Hzx1 = mdspan_t{&emdata.Hz[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny - 1  , Nz      }};
      const auto Hxy1 = mdspan_t{&emdata.Hx[0, Ny - PMLDepth, 0], std::extents{Nx      , PMLDepth, Nz - 1  }};
      const auto Hzy1 = mdspan_t{&emdata.Hz[0, Ny - PMLDepth, 0], std::extents{Nx - 1  , PMLDepth, Nz      }};
      const auto Hxz1 = mdspan_t{&emdata.Hx[0, 0, Nz - PMLDepth], std::extents{Nx      , Ny - 1  , PMLDepth}};
      const auto Hyz1 = mdspan_t{&emdata.Hy[0, 0, Nz - PMLDepth], std::extents{Nx - 1  , Ny      , PMLDepth}};

      const auto Eyx0 = mdspan_t{&emdata.Ey[0, 0, 0], std::extents{PMLDepth, Ny - 1  , Nz      }};
      const auto Ezx0 = mdspan_t{&emdata.Ez[0, 0, 0], std::extents{PMLDepth, Ny      , Nz - 1  }};
      const auto Exy0 = mdspan_t{&emdata.Ex[0, 0, 0], std::extents{Nx - 1  , PMLDepth, Nz      }};
      const auto Ezy0 = mdspan_t{&emdata.Ez[0, 0, 0], std::extents{Nx      , PMLDepth, Nz - 1  }};
      const auto Exz0 = mdspan_t{&emdata.Ex[0, 0, 0], std::extents{Nx - 1  , Ny      , PMLDepth}};
      const auto Eyz0 = mdspan_t{&emdata.Ey[0, 0, 0], std::extents{Nx      , Ny - 1  , PMLDepth}};

      const auto Eyx1 = mdspan_t{&emdata.Ey[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny - 1  , Nz      }};
      const auto Ezx1 = mdspan_t{&emdata.Ez[Nx - PMLDepth, 0, 0], std::extents{PMLDepth, Ny      , Nz - 1  }};
      const auto Exy1 = mdspan_t{&emdata.Ex[0, Ny - PMLDepth, 0], std::extents{Nx - 1  , PMLDepth, Nz      }};
      const auto Ezy1 = mdspan_t{&emdata.Ez[0, Ny - PMLDepth, 0], std::extents{Nx      , PMLDepth, Nz - 1  }};
      const auto Exz1 = mdspan_t{&emdata.Ex[0, 0, Nz - PMLDepth], std::extents{Nx - 1  , Ny      , PMLDepth}};
      const auto Eyz1 = mdspan_t{&emdata.Ey[0, 0, Nz - PMLDepth], std::extents{Nx      , Ny - 1  , PMLDepth}};

      constexpr auto Chxy = dt / (constants::mu0<double> * dz);
      constexpr auto Chxz = dt / (constants::mu0<double> * dy);
      constexpr auto Chyz = dt / (constants::mu0<double> * dx);

      X0Hy::apply(Hyx0, Ezx0, Chyz, emdata.x0_Hy);
      X1Hy::apply(Hyx1, Ezx1, Chyz, emdata.x1_Hy);
   
      X0Hz::apply(Hzx0, Eyx0, Chyz, emdata.x0_Hz);
      X1Hz::apply(Hzx1, Eyx1, Chyz, emdata.x1_Hz);
   
      Y0Hx::apply(Hxy0, Ezy0, Chxz, emdata.y0_Hx);
      Y1Hx::apply(Hxy1, Ezy1, Chxz, emdata.y1_Hx);
   
      Y0Hz::apply(Hzy0, Exy0, Chxz, emdata.y0_Hz);
      Y1Hz::apply(Hzy1, Exy1, Chxz, emdata.y1_Hz);
   
      Z0Hx::apply(Hxz0, Eyz0, Chxy, emdata.z0_Hx);
      Z1Hx::apply(Hxz1, Eyz1, Chxy, emdata.z1_Hx);
   
      Z0Hy::apply(Hyz0, Exz0, Chxy, emdata.z0_Hy);
      Z1Hy::apply(Hyz1, Exz1, Chxy, emdata.z1_Hy);
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
