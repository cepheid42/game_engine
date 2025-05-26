#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP

#include "program_params.hpp"
#include "em_data.hpp"
#include "bc_data.hpp"
#include "update_functors.hpp"
#include "bc_functors.hpp"
#include "diff_operators.hpp"

#include <algorithm>

namespace tf::electromagnetics
{
struct EMSolver
{
   using ex_func = FieldIntegrator<ExplicitUpdateFunctor<noop, backward_dz>>;
   using ey_func = FieldIntegrator<ExplicitUpdateFunctor<backward_dz, backward_dx>>;
   using ez_func = FieldIntegrator<ExplicitUpdateFunctor<backward_dx, noop>>;
   using hx_func = FieldIntegrator<ExplicitUpdateFunctor<forward_dz, noop>>;
   using hy_func = FieldIntegrator<ExplicitUpdateFunctor<forward_dx, forward_dz>>;
   using hz_func = FieldIntegrator<ExplicitUpdateFunctor<noop, forward_dx>>;

   // X-Faces
   using Ex_x0_bc = BCIntegrator<void>;
   using Ex_x1_bc = BCIntegrator<void>;
   using Ey_x0_bc = BCIntegrator<PMLFunctor<backward_dx, false, true>>;
   using Ey_x1_bc = BCIntegrator<PMLFunctor<backward_dx, true, true>>;
   using Ez_x0_bc = BCIntegrator<PMLFunctor<backward_dx, false, false>>;
   using Ez_x1_bc = BCIntegrator<PMLFunctor<backward_dx, true, false>>;

   using Hx_x0_bc = BCIntegrator<void>;
   using Hx_x1_bc = BCIntegrator<void>;
   using Hy_x0_bc = BCIntegrator<PMLFunctor<forward_dx, false, false>>;
   using Hy_x1_bc = BCIntegrator<PMLFunctor<forward_dx, true, false>>;
   using Hz_x0_bc = BCIntegrator<PMLFunctor<forward_dx, false, true>>;
   using Hz_x1_bc = BCIntegrator<PMLFunctor<forward_dx, true, true>>;

   // Y-Faces
   using Ex_y0_bc = BCIntegrator<void>;
   using Ex_y1_bc = BCIntegrator<void>;
   using Ey_y0_bc = BCIntegrator<void>;
   using Ey_y1_bc = BCIntegrator<void>;
   using Ez_y0_bc = BCIntegrator<void>;
   using Ez_y1_bc = BCIntegrator<void>;

   using Hx_y0_bc = BCIntegrator<void>;
   using Hx_y1_bc = BCIntegrator<void>;
   using Hy_y0_bc = BCIntegrator<void>;
   using Hy_y1_bc = BCIntegrator<void>;
   using Hz_y0_bc = BCIntegrator<void>;
   using Hz_y1_bc = BCIntegrator<void>;

   // Z-Faces
   using Ex_z0_bc = BCIntegrator<PMLFunctor<backward_dz, false, true>>;
   using Ex_z1_bc = BCIntegrator<PMLFunctor<backward_dz, true, true>>;
   using Ey_z0_bc = BCIntegrator<PMLFunctor<backward_dz, false, false>>;
   using Ey_z1_bc = BCIntegrator<PMLFunctor<backward_dz, true, false>>;
   using Ez_z0_bc = BCIntegrator<void>;
   using Ez_z1_bc = BCIntegrator<void>;

   using Hx_z0_bc = BCIntegrator<PMLFunctor<forward_dz, false, false>>;
   using Hx_z1_bc = BCIntegrator<PMLFunctor<forward_dz, true, false>>;
   using Hy_z0_bc = BCIntegrator<PMLFunctor<forward_dz, false, true>>;
   using Hy_z1_bc = BCIntegrator<PMLFunctor<forward_dz, true, true>>;
   using Hz_z0_bc = BCIntegrator<void>;
   using Hz_z1_bc = BCIntegrator<void>;


   explicit EMSolver() = delete;

   explicit EMSolver(const std::size_t nx, const std::size_t ny, const std::size_t nz, const compute_t cfl,
                     const compute_t   dt)
   : emdata(nx, ny, nz, cfl, dt),
     bcdata(this->emdata)
   {}

   void advance(const compute_t t) {
      updateH();
      updateHBCs();
      apply_srcs(t);
      updateE();
      updateEBCs();
      updateBhalf();   // for the particles and shit
      zero_currents(); // also for the particles, don't need last weeks currents
   }

   void updateE() {
      // todo: changed y-limits to {0, 0} otherwise y-loop never executes, Dy is now a NoOp
      ex_update(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexhz, emdata.Cexhy, emdata.Cjx,
                {0, 0, 0, 0, 1, 1});
      ey_update(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyhx, emdata.Ceyhz, emdata.Cjy,
                {1, 1, 0, 0, 1, 1});
      ez_update(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezhy, emdata.Cezhx, emdata.Cjz,
                {1, 1, 0, 0, 0, 0});
   }

   void updateH() {
      hx_update(emdata.Hx, emdata.Ey, emdata.Ez, emdata.empty, emdata.Chxh, emdata.Chxey, emdata.Chxez, emdata.empty,
                {0, 0, 0, 0, 0, 0});
      hy_update(emdata.Hy, emdata.Ez, emdata.Ex, emdata.empty, emdata.Chyh, emdata.Chyez, emdata.Chyex, emdata.empty,
                {0, 0, 0, 0, 0, 0});
      hz_update(emdata.Hz, emdata.Ex, emdata.Ey, emdata.empty, emdata.Chzh, emdata.Chzex, emdata.Chzey, emdata.empty,
                {0, 0, 0, 0, 0, 0});
   }

   void updateBhalf() {
      std::ranges::copy(emdata.Hx.begin(), emdata.Hx.end(), emdata.Bx.begin());
      std::ranges::copy(emdata.Hy.begin(), emdata.Hy.end(), emdata.By.begin());
      std::ranges::copy(emdata.Hz.begin(), emdata.Hz.end(), emdata.Bz.begin());

      hx_update(emdata.Bx, emdata.Ey, emdata.Ez, emdata.empty, emdata.Chxh, emdata.Chxey2, emdata.Chxez2, emdata.empty,
                {0, 0, 0, 0, 0, 0});
      hy_update(emdata.By, emdata.Ez, emdata.Ex, emdata.empty, emdata.Chyh, emdata.Chyez2, emdata.Chyex2, emdata.empty,
                {0, 0, 0, 0, 0, 0});
      hz_update(emdata.Bz, emdata.Ex, emdata.Ey, emdata.empty, emdata.Chzh, emdata.Chzex2, emdata.Chzey2, emdata.empty,
                {0, 0, 0, 0, 0, 0});

      // std::ranges::for_each(emdata.Bx.begin(), emdata.Bx.end(), [](auto& x) { x *= constants::mu0<compute_t>;});
      // std::ranges::for_each(emdata.By.begin(), emdata.By.end(), [](auto& x) { x *= constants::mu0<compute_t>; });
      // std::ranges::for_each(emdata.Bz.begin(), emdata.Bz.end(), [](auto& x) { x *= constants::mu0<compute_t>; });

      for (std::size_t i = 0; i < emdata.Bx.dims()[0]; i++) {
         for (std::size_t j = 0; j < emdata.Bx.dims()[1]; j++) {
            for (std::size_t k = 0; k < emdata.Bx.dims()[2]; k++) {
               emdata.Bx(i, j, k) = emdata.Bx(i, j, k) * constants::mu0<compute_t> + emdata.Bx_app(i, j, k);
            }
         }
      }

      for (std::size_t i = 0; i < emdata.By.dims()[0]; i++) {
         for (std::size_t j = 0; j < emdata.By.dims()[1]; j++) {
            for (std::size_t k = 0; k < emdata.By.dims()[2]; k++) {
               emdata.By(i, j, k) = emdata.By(i, j, k) * constants::mu0<compute_t> + emdata.By_app(i, j, k);
            }
         }
      }

      for (std::size_t i = 0; i < emdata.Bz.dims()[0]; i++) {
         for (std::size_t j = 0; j < emdata.Bz.dims()[1]; j++) {
            for (std::size_t k = 0; k < emdata.Bz.dims()[2]; k++) {
               emdata.Bz(i, j, k) = emdata.Bz(i, j, k) * constants::mu0<compute_t> + emdata.Bz_app(i, j, k);
            }
         }
      }
   }

   void updateEBCs() {
      Ey_x0(emdata.Ey, emdata.Hz, emdata.Ceyhz, bcdata.x0.Ey);
      Ey_x1(emdata.Ey, emdata.Hz, emdata.Ceyhz, bcdata.x1.Ey);
      Ez_x0(emdata.Ez, emdata.Hy, emdata.Cezhy, bcdata.x0.Ez);
      Ez_x1(emdata.Ez, emdata.Hy, emdata.Cezhy, bcdata.x1.Ez);

      Ex_y0(emdata.Ex, emdata.Hz, emdata.Cexhz, bcdata.y0.Ex);
      Ex_y1(emdata.Ex, emdata.Hz, emdata.Cexhz, bcdata.y1.Ex);
      Ez_y0(emdata.Ez, emdata.Hx, emdata.Cezhx, bcdata.y0.Ez);
      Ez_y1(emdata.Ez, emdata.Hx, emdata.Cezhx, bcdata.y1.Ez);

      Ex_z0(emdata.Ex, emdata.Hy, emdata.Cexhy, bcdata.z0.Ex);
      Ex_z1(emdata.Ex, emdata.Hy, emdata.Cexhy, bcdata.z1.Ex);
      Ey_z0(emdata.Ey, emdata.Hx, emdata.Ceyhx, bcdata.z0.Ey);
      Ey_z1(emdata.Ey, emdata.Hx, emdata.Ceyhx, bcdata.z1.Ey);
   }

   void updateHBCs() {
      Hy_x0(emdata.Hy, emdata.Ez, emdata.Chyez, bcdata.x0.Hy);
      Hy_x1(emdata.Hy, emdata.Ez, emdata.Chyez, bcdata.x1.Hy);
      Hz_x0(emdata.Hz, emdata.Ey, emdata.Chzey, bcdata.x0.Hz);
      Hz_x1(emdata.Hz, emdata.Ey, emdata.Chzey, bcdata.x1.Hz);

      Hx_y0(emdata.Hx, emdata.Ez, emdata.Chxez, bcdata.y0.Hx);
      Hx_y1(emdata.Hx, emdata.Ez, emdata.Chxez, bcdata.y1.Hx);
      Hz_y0(emdata.Hz, emdata.Ex, emdata.Chzex, bcdata.y0.Hz);
      Hz_y1(emdata.Hz, emdata.Ex, emdata.Chzex, bcdata.y1.Hz);

      Hx_z0(emdata.Hx, emdata.Ey, emdata.Chxey, bcdata.z0.Hx);
      Hx_z1(emdata.Hx, emdata.Ey, emdata.Chxey, bcdata.z1.Hx);
      Hy_z0(emdata.Hy, emdata.Ex, emdata.Chyex, bcdata.z0.Hy);
      Hy_z1(emdata.Hy, emdata.Ex, emdata.Chyex, bcdata.z1.Hy);
   }

   void apply_srcs(const compute_t t) const {
      for (const auto& src: emdata.srcs)
      {
         src.apply(t);
      }
   }

   void zero_currents() {
      std::ranges::fill(emdata.Jx, 0.0_fp);
      std::ranges::fill(emdata.Jy, 0.0_fp);
      std::ranges::fill(emdata.Jz, 0.0_fp);
   }

   EMData emdata;
   BCData bcdata;

   ex_func ex_update{};
   ey_func ey_update{};
   ez_func ez_update{};
   hx_func hx_update{};
   hy_func hy_update{};
   hz_func hz_update{};

   Ey_x0_bc Ey_x0{};
   Ey_x1_bc Ey_x1{};
   Ez_x0_bc Ez_x0{};
   Ez_x1_bc Ez_x1{};

   Hy_x0_bc Hy_x0{};
   Hy_x1_bc Hy_x1{};
   Hz_x0_bc Hz_x0{};
   Hz_x1_bc Hz_x1{};

   Ex_y0_bc Ex_y0{};
   Ex_y1_bc Ex_y1{};
   Ez_y0_bc Ez_y0{};
   Ez_y1_bc Ez_y1{};

   Hx_y0_bc Hx_y0{};
   Hx_y1_bc Hx_y1{};
   Hz_y0_bc Hz_y0{};
   Hz_y1_bc Hz_y1{};

   Ex_z0_bc Ex_z0{};
   Ex_z1_bc Ex_z1{};
   Ey_z0_bc Ey_z0{};
   Ey_z1_bc Ey_z1{};

   Hx_z0_bc Hx_z0{};
   Hx_z1_bc Hx_z1{};
   Hy_z0_bc Hy_z0{};
   Hy_z1_bc Hy_z1{};
};
} // end namespace tf::electromagnetics

#endif //EM_SOLVER_HPP
