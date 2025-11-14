#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "em_definitions.hpp"
#include "em_data.hpp"

namespace tf::electromagnetics {

template<typename UpdateFunc>
struct FieldIntegrator {
   using offset_t = std::array<std::size_t, 6>;

   static void apply(auto& f, const auto& d1, const auto& d2, const auto& src,
                     const auto& c_f, const auto& c_d1, const auto& c_d2, const auto& c_src,
                     const offset_t& offsets)
   {
      #pragma omp parallel for simd collapse(3) num_threads(nThreads)
      for (std::size_t i = offsets[0]; i < f.nx() - offsets[1]; ++i) {
         for (std::size_t j = offsets[2]; j < f.ny() - offsets[3]; ++j) {
            for (std::size_t k = offsets[4]; k < f.nz() - offsets[5]; ++k) {
               UpdateFunc::apply(f, d1, d2, src, c_f, c_d1, c_d2, c_src, i, j, k);
            } // end for k
         } // end for j
      } // end for i
   } // end operator()
}; // end struct FieldIntegrator

struct EMSolver {
   static void advance(auto& emdata, const auto t) requires(em_enabled) {
      updateH();
      // updateHBCs();
      // updateJBCs();
      // apply_srcs(t);
      updateE();
      // updateEBCs();
      // particle_correction(); // for the particles and shit
      // zero_currents();       // also for the particles, don't need last week's currents
   }

   static void advance(auto&, const auto) requires (!em_enabled) {}

   void updateE(auto& emdata) {
      ExUpdate::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Ceyz, emdata.Cexy, emdata.Cj, Ex_offsets);
      EyUpdate::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Cexy, emdata.Ceyz, emdata.Cj, Ey_offsets);
      EzUpdate::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceyz, emdata.Cexz, emdata.Cj, Ez_offsets);
   }

   void updateH(auto& emdata) {
      HxUpdate::apply(emdata.Hx, emdata.Ey, emdata.Ez, emdata.empty, emdata.Chxy, emdata.Chxz, emdata.empty, {0, 0, 0, 0, 0, 0});
      HyUpdate::apply(emdata.Hy, emdata.Ez, emdata.Ex, emdata.empty, emdata.Chyz, emdata.Chxy, emdata.empty, {0, 0, 0, 0, 0, 0});
      HzUpdate::apply(emdata.Hz, emdata.Ex, emdata.Ey, emdata.empty, emdata.Chxz, emdata.Chyz, emdata.empty, {0, 0, 0, 0, 0, 0});
   }

   // void particle_correction() {
   //    std::ranges::copy(emdata.Hx, emdata.Bx.begin());
   //    std::ranges::copy(emdata.Hy, emdata.By.begin());
   //    std::ranges::copy(emdata.Hz, emdata.Bz.begin());
   //
   //    HxUpdate::apply(emdata.Bx, emdata.Ey, emdata.Ez, emdata.empty, emdata.Chxy2, emdata.Chxz2, emdata.empty, {0, 0, 0, 0, 0, 0});
   //    HyUpdate::apply(emdata.By, emdata.Ez, emdata.Ex, emdata.empty, emdata.Chyz2, emdata.Chxy2, emdata.empty, {0, 0, 0, 0, 0, 0});
   //    HzUpdate::apply(emdata.Bz, emdata.Ex, emdata.Ey, emdata.empty, emdata.Chxz2, emdata.Chyz2, emdata.empty, {0, 0, 0, 0, 0, 0});
   //
   //    #pragma omp parallel num_threads(nThreads)
   //    {
   //       #pragma omp for simd
   //       for (std::size_t i = 0; i < emdata.Ex.size(); i++) {
   //          emdata.Ex_total[i] = emdata.Ex[i] + emdata.Ex_app[i];
   //       }
   //
   //       #pragma omp for simd
   //       for (std::size_t i = 0; i < emdata.Ey.size(); i++) {
   //          emdata.Ey_total[i] = emdata.Ey[i] + emdata.Ey_app[i];
   //       }
   //
   //       #pragma omp for simd
   //       for (std::size_t i = 0; i < emdata.Ez.size(); i++) {
   //          emdata.Ez_total[i] = emdata.Ez[i] + emdata.Ez_app[i];
   //       }
   //
   //       #pragma omp for simd
   //       for (std::size_t i = 0; i < emdata.Bx.size(); i++) {
   //          emdata.Bx_total[i] = emdata.Bx[i] * constants::mu0<double> + emdata.Bx_app[i];
   //       }
   //
   //       #pragma omp for simd
   //       for (std::size_t i = 0; i < emdata.By.size(); i++) {
   //          emdata.By_total[i] = emdata.By[i] * constants::mu0<double> + emdata.By_app[i];
   //       }
   //
   //       #pragma omp for simd
   //       for (std::size_t i = 0; i < emdata.Bz.size(); i++) {
   //          emdata.Bz_total[i] = emdata.Bz[i] * constants::mu0<double> + emdata.Bz_app[i];
   //       }
   //    } // end omp parallel
   // }
   //
   // void updateEBCs() {
   //    X0BC::Ey::apply(emdata.Ey, emdata.Hz, emdata.Ceyhz, bcdata.x0.Ey);
   //    X1BC::Ey::apply(emdata.Ey, emdata.Hz, emdata.Ceyhz, bcdata.x1.Ey);
   //
   //    X0BC::Ez::apply(emdata.Ez, emdata.Hy, emdata.Cezhy, bcdata.x0.Ez);
   //    X1BC::Ez::apply(emdata.Ez, emdata.Hy, emdata.Cezhy, bcdata.x1.Ez);
   //
   //    Y0BC::Ex::apply(emdata.Ex, emdata.Hz, emdata.Cexhz, bcdata.y0.Ex);
   //    Y1BC::Ex::apply(emdata.Ex, emdata.Hz, emdata.Cexhz, bcdata.y1.Ex);
   //
   //    Y0BC::Ez::apply(emdata.Ez, emdata.Hx, emdata.Cezhx, bcdata.y0.Ez);
   //    Y1BC::Ez::apply(emdata.Ez, emdata.Hx, emdata.Cezhx, bcdata.y1.Ez);
   //
   //    Z0BC::Ex::apply(emdata.Ex, emdata.Hy, emdata.Cexhy, bcdata.z0.Ex);
   //    Z1BC::Ex::apply(emdata.Ex, emdata.Hy, emdata.Cexhy, bcdata.z1.Ex);
   //
   //    Z0BC::Ey::apply(emdata.Ey, emdata.Hx, emdata.Ceyhx, bcdata.z0.Ey);
   //    Z1BC::Ey::apply(emdata.Ey, emdata.Hx, emdata.Ceyhx, bcdata.z1.Ey);
   // }
   //
   // void updateHBCs() {
   //    X0BC::Hy::apply(emdata.Hy, emdata.Ez, emdata.Chyez, bcdata.x0.Hy);
   //    X1BC::Hy::apply(emdata.Hy, emdata.Ez, emdata.Chyez, bcdata.x1.Hy);
   //
   //    X0BC::Hz::apply(emdata.Hz, emdata.Ey, emdata.Chzey, bcdata.x0.Hz);
   //    X1BC::Hz::apply(emdata.Hz, emdata.Ey, emdata.Chzey, bcdata.x1.Hz);
   //
   //    Y0BC::Hx::apply(emdata.Hx, emdata.Ez, emdata.Chxez, bcdata.y0.Hx);
   //    Y1BC::Hx::apply(emdata.Hx, emdata.Ez, emdata.Chxez, bcdata.y1.Hx);
   //
   //    Y0BC::Hz::apply(emdata.Hz, emdata.Ex, emdata.Chzex, bcdata.y0.Hz);
   //    Y1BC::Hz::apply(emdata.Hz, emdata.Ex, emdata.Chzex, bcdata.y1.Hz);
   //
   //    Z0BC::Hx::apply(emdata.Hx, emdata.Ey, emdata.Chxey, bcdata.z0.Hx);
   //    Z1BC::Hx::apply(emdata.Hx, emdata.Ey, emdata.Chxey, bcdata.z1.Hx);
   //
   //    Z0BC::Hy::apply(emdata.Hy, emdata.Ex, emdata.Chyex, bcdata.z0.Hy);
   //    Z1BC::Hy::apply(emdata.Hy, emdata.Ex, emdata.Chyex, bcdata.z1.Hy);
   // }
   //
   // void updateJBCs() {
   //    // Only used for periodic BCs
   //    X0BC::Jy::apply(emdata.Jy, emdata.empty, emdata.empty, bcdata.x0.Jy);
   //    X0BC::Jz::apply(emdata.Jz, emdata.empty, emdata.empty, bcdata.x0.Jz);
   //
   //    Y0BC::Jx::apply(emdata.Jx, emdata.empty, emdata.empty, bcdata.y0.Jx);
   //    Y0BC::Jz::apply(emdata.Jz, emdata.empty, emdata.empty, bcdata.y0.Jz);
   //
   //    Z0BC::Jx::apply(emdata.Jx, emdata.empty, emdata.empty, bcdata.z0.Jx);
   //    Z0BC::Jy::apply(emdata.Jy, emdata.empty, emdata.empty, bcdata.z0.Jy);
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
