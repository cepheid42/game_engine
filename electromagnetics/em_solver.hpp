#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "em_definitions.hpp"
#include "em_data.hpp"

namespace tf::electromagnetics {


struct EMSolver {
   using offset_t = std::array<std::size_t, 6>;

   using ExUpdate = TypeListAt<0, field_t>;
   using EyUpdate = TypeListAt<1, field_t>;
   using EzUpdate = TypeListAt<2, field_t>;
   using HxUpdate = TypeListAt<3, field_t>;
   using HyUpdate = TypeListAt<4, field_t>;
   using HzUpdate = TypeListAt<5, field_t>;

   using X0BC = TypeListAt<0, boundary_t>;
   using X1BC = TypeListAt<1, boundary_t>;
   using Y0BC = TypeListAt<2, boundary_t>;
   using Y1BC = TypeListAt<3, boundary_t>;
   using Z0BC = TypeListAt<4, boundary_t>;
   using Z1BC = TypeListAt<5, boundary_t>;

   // todo: support other 2D modes and 1D
   static constexpr offset_t Ex_offsets = {0,            0,            !y_collapsed, !y_collapsed, !z_collapsed, !z_collapsed};
   static constexpr offset_t Ey_offsets = {!x_collapsed, !x_collapsed, 0,            0,            !z_collapsed, !z_collapsed};
   static constexpr offset_t Ez_offsets = {!x_collapsed, !x_collapsed, !y_collapsed, !y_collapsed, 0,            0};

   explicit EMSolver() = delete;

   explicit EMSolver(const std::size_t nx, const std::size_t ny, const std::size_t nz)
   : emdata(nx, ny, nz),
     bcdata(this->emdata)
   {}

   void advance(const auto t) requires(em_enabled) {
      for (auto i = 0zu; i < em_subcycles; i++) {
         updateH();
         updateHBCs();
         updateJBCs();
         apply_srcs(t);
         updateE();
         updateEBCs();
      }

      if constexpr (push_enabled) {
         particle_correction(); // for the particles and shit
      }

      if constexpr (jdep_enabled) {
         zero_currents();       // also for the particles, don't need last week's currents
      }
   }

   static void advance(const auto) requires (!em_enabled) {}

   void updateE() {
      ExUpdate::apply(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexhz, emdata.Cexhy, emdata.Cjx, Ex_offsets);
      EyUpdate::apply(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyhx, emdata.Ceyhz, emdata.Cjy, Ey_offsets);
      EzUpdate::apply(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezhy, emdata.Cezhx, emdata.Cjz, Ez_offsets);
   }

   void updateH() {
      HxUpdate::apply(emdata.Hx, emdata.Ey, emdata.Ez, emdata.empty, emdata.Chxh, emdata.Chxey, emdata.Chxez, emdata.empty, {0, 0, 0, 0, 0, 0});
      HyUpdate::apply(emdata.Hy, emdata.Ez, emdata.Ex, emdata.empty, emdata.Chyh, emdata.Chyez, emdata.Chyex, emdata.empty, {0, 0, 0, 0, 0, 0});
      HzUpdate::apply(emdata.Hz, emdata.Ex, emdata.Ey, emdata.empty, emdata.Chzh, emdata.Chzex, emdata.Chzey, emdata.empty, {0, 0, 0, 0, 0, 0});
   }

   void particle_correction() {
      std::ranges::copy(emdata.Hx, emdata.Bx.begin());
      std::ranges::copy(emdata.Hy, emdata.By.begin());
      std::ranges::copy(emdata.Hz, emdata.Bz.begin());

      HxUpdate::apply(emdata.Bx, emdata.Ey, emdata.Ez, emdata.empty, emdata.Chxh, emdata.Chxey2, emdata.Chxez2, emdata.empty, {0, 0, 0, 0, 0, 0});
      HyUpdate::apply(emdata.By, emdata.Ez, emdata.Ex, emdata.empty, emdata.Chyh, emdata.Chyez2, emdata.Chyex2, emdata.empty, {0, 0, 0, 0, 0, 0});
      HzUpdate::apply(emdata.Bz, emdata.Ex, emdata.Ey, emdata.empty, emdata.Chzh, emdata.Chzex2, emdata.Chzey2, emdata.empty, {0, 0, 0, 0, 0, 0});

      #pragma omp parallel num_threads(nThreads)
      {
         #pragma omp for
         for (std::size_t i = 0; i < emdata.Ex.size(); i++) {
            emdata.Ex_total[i] = emdata.Ex[i] + emdata.Ex_app[i];
         }

         #pragma omp for
         for (std::size_t i = 0; i < emdata.Ey.size(); i++) {
            emdata.Ey_total[i] = emdata.Ey[i] + emdata.Ey_app[i];
         }

         #pragma omp for
         for (std::size_t i = 0; i < emdata.Ez.size(); i++) {
            emdata.Ez_total[i] = emdata.Ez[i] + emdata.Ez_app[i];
         }

         #pragma omp for
         for (std::size_t i = 0; i < emdata.Bx.size(); i++) {
            emdata.Bx_total[i] = emdata.Bx[i] * constants::mu0 + emdata.Bx_app[i];
         }

         #pragma omp for
         for (std::size_t i = 0; i < emdata.By.size(); i++) {
            emdata.By_total[i] = emdata.By[i] * constants::mu0 + emdata.By_app[i];
         }

         #pragma omp for
         for (std::size_t i = 0; i < emdata.Bz.size(); i++) {
            emdata.Bz_total[i] = emdata.Bz[i] * constants::mu0 + emdata.Bz_app[i];
         }
      } // end omp parallel
   }

   void updateEBCs() {
      X0BC::Ey::apply(emdata.Ey, emdata.Hz, emdata.Ceyhz, bcdata.x0.Ey);
      X1BC::Ey::apply(emdata.Ey, emdata.Hz, emdata.Ceyhz, bcdata.x1.Ey);

      X0BC::Ez::apply(emdata.Ez, emdata.Hy, emdata.Cezhy, bcdata.x0.Ez);
      X1BC::Ez::apply(emdata.Ez, emdata.Hy, emdata.Cezhy, bcdata.x1.Ez);

      Y0BC::Ex::apply(emdata.Ex, emdata.Hz, emdata.Cexhz, bcdata.y0.Ex);
      Y1BC::Ex::apply(emdata.Ex, emdata.Hz, emdata.Cexhz, bcdata.y1.Ex);

      Y0BC::Ez::apply(emdata.Ez, emdata.Hx, emdata.Cezhx, bcdata.y0.Ez);
      Y1BC::Ez::apply(emdata.Ez, emdata.Hx, emdata.Cezhx, bcdata.y1.Ez);

      Z0BC::Ex::apply(emdata.Ex, emdata.Hy, emdata.Cexhy, bcdata.z0.Ex);
      Z1BC::Ex::apply(emdata.Ex, emdata.Hy, emdata.Cexhy, bcdata.z1.Ex);

      Z0BC::Ey::apply(emdata.Ey, emdata.Hx, emdata.Ceyhx, bcdata.z0.Ey);
      Z1BC::Ey::apply(emdata.Ey, emdata.Hx, emdata.Ceyhx, bcdata.z1.Ey);
   }

   void updateHBCs() {
      X0BC::Hy::apply(emdata.Hy, emdata.Ez, emdata.Chyez, bcdata.x0.Hy);
      X1BC::Hy::apply(emdata.Hy, emdata.Ez, emdata.Chyez, bcdata.x1.Hy);

      X0BC::Hz::apply(emdata.Hz, emdata.Ey, emdata.Chzey, bcdata.x0.Hz);
      X1BC::Hz::apply(emdata.Hz, emdata.Ey, emdata.Chzey, bcdata.x1.Hz);

      Y0BC::Hx::apply(emdata.Hx, emdata.Ez, emdata.Chxez, bcdata.y0.Hx);
      Y1BC::Hx::apply(emdata.Hx, emdata.Ez, emdata.Chxez, bcdata.y1.Hx);

      Y0BC::Hz::apply(emdata.Hz, emdata.Ex, emdata.Chzex, bcdata.y0.Hz);
      Y1BC::Hz::apply(emdata.Hz, emdata.Ex, emdata.Chzex, bcdata.y1.Hz);

      Z0BC::Hx::apply(emdata.Hx, emdata.Ey, emdata.Chxey, bcdata.z0.Hx);
      Z1BC::Hx::apply(emdata.Hx, emdata.Ey, emdata.Chxey, bcdata.z1.Hx);

      Z0BC::Hy::apply(emdata.Hy, emdata.Ex, emdata.Chyex, bcdata.z0.Hy);
      Z1BC::Hy::apply(emdata.Hy, emdata.Ex, emdata.Chyex, bcdata.z1.Hy);
   }

   void updateJBCs() {
      // Only used for periodic BCs
      X0BC::Jy::apply(emdata.Jy, emdata.empty, emdata.empty, bcdata.x0.Jy);
      X0BC::Jz::apply(emdata.Jz, emdata.empty, emdata.empty, bcdata.x0.Jz);

      Y0BC::Jx::apply(emdata.Jx, emdata.empty, emdata.empty, bcdata.y0.Jx);
      Y0BC::Jz::apply(emdata.Jz, emdata.empty, emdata.empty, bcdata.y0.Jz);

      Z0BC::Jx::apply(emdata.Jx, emdata.empty, emdata.empty, bcdata.z0.Jx);
      Z0BC::Jy::apply(emdata.Jy, emdata.empty, emdata.empty, bcdata.z0.Jy);
   }


   void apply_srcs(const double t) {
      for (const auto& src: emdata.srcs) {
         src.apply(t);
      }

      for (const auto& src: emdata.beams) {
         src.apply(t);
      }
   }

   void zero_currents() {
      std::ranges::fill(emdata.Jx, 0.0);
      std::ranges::fill(emdata.Jy, 0.0);
      std::ranges::fill(emdata.Jz, 0.0);
   }

   emdata_t emdata;
   bcdata_t bcdata;
}; // end struct EMSolver

using emsolver_t = EMSolver;
} // end namespace tf::electromagnetics

#endif //EM_SOLVER_HPP
