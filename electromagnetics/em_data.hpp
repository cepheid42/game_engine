#ifndef EM_DATA_HPP
#define EM_DATA_HPP

#include "program_params.hpp"
#include "sources.hpp"
#include "array.hpp"

namespace tf::electromagnetics {
   struct EMData {
      EMData() = delete;

      EMData(const std::size_t nx, const std::size_t ny, const std::size_t nz, const compute_t cfl, const compute_t dt)
      : Ex(nx - 1, ny, nz),
        Jx(nx - 1, ny, nz),
        Cexe(nx - 1, ny, nz),
        Cexhy(nx - 1, ny, nz),
        Cexhz(nx - 1, ny, nz),
        Cjx(nx - 1, ny, nz),
        Ey(nx, ny - 1, nz),
        Jy(nx, ny - 1, nz),
        Ceye(nx, ny - 1, nz),
        Ceyhx(nx, ny - 1, nz),
        Ceyhz(nx, ny - 1, nz),
        Cjy(nx, ny - 1, nz),
        Ez(nx, ny, nz - 1),
        Jz(nx, ny, nz - 1),
        Ceze(nx, ny, nz - 1),
        Cezhx(nx, ny, nz - 1),
        Cezhy(nx, ny, nz - 1),
        Cjz(nx, ny, nz - 1),
        Hx(nx, ny - 1, nz - 1),
        Chxh(nx, ny - 1, nz - 1),
        Chxey(nx, ny - 1, nz - 1),
        Chxez(nx, ny - 1, nz - 1),
        Hy(nx - 1, ny, nz - 1),
        Chyh(nx - 1, ny, nz - 1),
        Chyex(nx - 1, ny, nz - 1),
        Chyez(nx - 1, ny, nz - 1),
        Hz(nx - 1, ny - 1, nz),
        Chzh(nx - 1, ny - 1, nz),
        Chzex(nx - 1, ny - 1, nz),
        Chzey(nx - 1, ny - 1, nz),
        Bx(nx, ny - 1, nz - 1),
        Chxey2(nx, ny - 1, nz - 1),
        Chxez2(nx, ny - 1, nz - 1),
        By(nx - 1, ny, nz - 1),
        Chyex2(nx - 1, ny, nz - 1),
        Chyez2(nx - 1, ny, nz - 1),
        Bz(nx - 1, ny - 1, nz),
        Chzex2(nx - 1, ny - 1, nz),
        Chzey2(nx - 1, ny - 1, nz),
        Ex_app(nx - 1, ny, nz),
        Ey_app(nx, ny - 1, nz),
        Ez_app(nx, ny, nz - 1),
        Bx_app(nx, ny - 1, nz - 1),
        By_app(nx - 1, ny, nz - 1),
        Bz_app(nx - 1, ny - 1, nz),
        Ex_total(nx - 1, ny, nz),
        Ey_total(nx, ny - 1, nz),
        Ez_total(nx, ny, nz - 1),
        Bx_total(nx, ny - 1, nz - 1),
        By_total(nx - 1, ny, nz - 1),
        Bz_total(nx - 1, ny - 1, nz)
      {
         init_coefficients(cfl, dt);
      }

      void init_coefficients(const compute_t cfl, const compute_t dt) {
         const auto e_coeff = cfl * constants::eta0<compute_t>;
         const auto h_coeff = cfl / constants::eta0<compute_t>;
         const auto j_coeff = dt / constants::eps0<compute_t>;

         // todo: add loss terms
         Cexe.fill(1.0_fp);
         Ceye.fill(1.0_fp);
         Ceze.fill(1.0_fp);
         Chxh.fill(1.0_fp);
         Chyh.fill(1.0_fp);
         Chzh.fill(1.0_fp);

         Cexhy.fill(e_coeff);
         Cexhz.fill(e_coeff);
         Ceyhx.fill(e_coeff);
         Ceyhz.fill(e_coeff);
         Cezhx.fill(e_coeff);
         Cezhy.fill(e_coeff);

         Cjx.fill(j_coeff);
         Cjy.fill(j_coeff);
         Cjz.fill(j_coeff);

         Chxey.fill(h_coeff);
         Chxez.fill(h_coeff);
         Chyex.fill(h_coeff);
         Chyez.fill(h_coeff);
         Chzex.fill(h_coeff);
         Chzey.fill(h_coeff);

         Chxey2.fill(h_coeff / 2.0);
         Chxez2.fill(h_coeff / 2.0);
         Chyex2.fill(h_coeff / 2.0);
         Chyez2.fill(h_coeff / 2.0);
         Chzex2.fill(h_coeff / 2.0);
         Chzey2.fill(h_coeff / 2.0);
      }

      Array3D<compute_t> Ex;
      Array3D<compute_t> Jx;
      Array3D<compute_t> Cexe;
      Array3D<compute_t> Cexhy;
      Array3D<compute_t> Cexhz;
      Array3D<compute_t> Cjx;

      Array3D<compute_t> Ey;
      Array3D<compute_t> Jy;
      Array3D<compute_t> Ceye;
      Array3D<compute_t> Ceyhx;
      Array3D<compute_t> Ceyhz;
      Array3D<compute_t> Cjy;

      Array3D<compute_t> Ez;
      Array3D<compute_t> Jz;
      Array3D<compute_t> Ceze;
      Array3D<compute_t> Cezhx;
      Array3D<compute_t> Cezhy;
      Array3D<compute_t> Cjz;

      Array3D<compute_t> Hx;
      Array3D<compute_t> Chxh;
      Array3D<compute_t> Chxey;
      Array3D<compute_t> Chxez;

      Array3D<compute_t> Hy;
      Array3D<compute_t> Chyh;
      Array3D<compute_t> Chyex;
      Array3D<compute_t> Chyez;

      Array3D<compute_t> Hz;
      Array3D<compute_t> Chzh;
      Array3D<compute_t> Chzex;
      Array3D<compute_t> Chzey;

      Array3D<compute_t> Bx;
      Array3D<compute_t> Chxey2;
      Array3D<compute_t> Chxez2;

      Array3D<compute_t> By;
      Array3D<compute_t> Chyex2;
      Array3D<compute_t> Chyez2;

      Array3D<compute_t> Bz;
      Array3D<compute_t> Chzex2;
      Array3D<compute_t> Chzey2;

      Array3D<compute_t> Ex_app;
      Array3D<compute_t> Ey_app;
      Array3D<compute_t> Ez_app;

      Array3D<compute_t> Bx_app;
      Array3D<compute_t> By_app;
      Array3D<compute_t> Bz_app;

      Array3D<compute_t> Ex_total;
      Array3D<compute_t> Ey_total;
      Array3D<compute_t> Ez_total;

      Array3D<compute_t> Bx_total;
      Array3D<compute_t> By_total;
      Array3D<compute_t> Bz_total;

      std::vector<GaussianBeam> srcs{};

      Array3D<void> empty{}; // for the shits and possibly also some giggles...
   };
} // end namespace tf::electromagnetics


#endif //EM_DATA_HPP
