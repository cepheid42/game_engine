#ifndef EM_DATA_HPP
#define EM_DATA_HPP

#include "array.hpp"
#include "program_params.hpp"
#include "sources.hpp"

#include <adios2.h>
#include <unordered_map>

namespace tf::electromagnetics {
struct EMData {
   EMData() = delete;

   EMData(const std::size_t nx, const std::size_t ny, const std::size_t nz)
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
      init_coefficients();
      load_applied_fields();
   }

   void init_coefficients() {
      constexpr auto e_coeff = dt / constants::eps0;
      constexpr auto h_coeff = dt / constants::mu0;

      // todo: add loss terms
      Cexe.fill(1.0);
      Ceye.fill(1.0);
      Ceze.fill(1.0);
      Chxh.fill(1.0);
      Chyh.fill(1.0);
      Chzh.fill(1.0);

      Cexhy.fill(e_coeff / dz);
      Cexhz.fill(e_coeff / dy);
      Ceyhx.fill(e_coeff / dz);
      Ceyhz.fill(e_coeff / dx);
      Cezhx.fill(e_coeff / dy);
      Cezhy.fill(e_coeff / dx);

      Cjx.fill(e_coeff);
      Cjy.fill(e_coeff);
      Cjz.fill(e_coeff);

      Chxey.fill(h_coeff / dz);
      Chxez.fill(h_coeff / dy);
      Chyex.fill(h_coeff / dz);
      Chyez.fill(h_coeff / dx);
      Chzex.fill(h_coeff / dy);
      Chzey.fill(h_coeff / dx);

      Chxey2.fill(0.5 * h_coeff / dz);
      Chxez2.fill(0.5 * h_coeff / dy);
      Chyex2.fill(0.5 * h_coeff / dz);
      Chyez2.fill(0.5 * h_coeff / dx);
      Chzex2.fill(0.5 * h_coeff / dy);
      Chzey2.fill(0.5 * h_coeff / dx);
   }

   void load_applied_fields() {
      adios2::ADIOS adios{};
      adios2::IO io = adios.DeclareIO("EMFieldsExternal");
      adios2::Engine reader = io.Open(applied_fields_path, adios2::Mode::Read);

      reader.BeginStep();

      const auto bpEx = io.InquireVariable<double>("Ex");
      const auto bpEy = io.InquireVariable<double>("Ey");
      const auto bpEz = io.InquireVariable<double>("Ez");
      const auto bpBx = io.InquireVariable<double>("Bx");
      const auto bpBy = io.InquireVariable<double>("By");
      const auto bpBz = io.InquireVariable<double>("Bz");

      if (bpEx) { reader.Get(bpEx, Ex_app.data(), adios2::Mode::Sync); }
      if (bpEy) { reader.Get(bpEy, Ey_app.data(), adios2::Mode::Sync); }
      if (bpEy) { reader.Get(bpEz, Ez_app.data(), adios2::Mode::Sync); }
      if (bpBx) { reader.Get(bpBx, Bx_app.data(), adios2::Mode::Sync); }
      if (bpBy) { reader.Get(bpBy, By_app.data(), adios2::Mode::Sync); }
      if (bpBz) { reader.Get(bpBz, Bz_app.data(), adios2::Mode::Sync); }

      reader.EndStep();
      reader.Close();
   }

   Array3D<double> Ex;
   Array3D<double> Jx;
   Array3D<double> Cexe;
   Array3D<double> Cexhy;
   Array3D<double> Cexhz;
   Array3D<double> Cjx;

   Array3D<double> Ey;
   Array3D<double> Jy;
   Array3D<double> Ceye;
   Array3D<double> Ceyhx;
   Array3D<double> Ceyhz;
   Array3D<double> Cjy;

   Array3D<double> Ez;
   Array3D<double> Jz;
   Array3D<double> Ceze;
   Array3D<double> Cezhx;
   Array3D<double> Cezhy;
   Array3D<double> Cjz;

   Array3D<double> Hx;
   Array3D<double> Chxh;
   Array3D<double> Chxey;
   Array3D<double> Chxez;

   Array3D<double> Hy;
   Array3D<double> Chyh;
   Array3D<double> Chyex;
   Array3D<double> Chyez;

   Array3D<double> Hz;
   Array3D<double> Chzh;
   Array3D<double> Chzex;
   Array3D<double> Chzey;

   Array3D<double> Bx;
   Array3D<double> Chxey2;
   Array3D<double> Chxez2;

   Array3D<double> By;
   Array3D<double> Chyex2;
   Array3D<double> Chyez2;

   Array3D<double> Bz;
   Array3D<double> Chzex2;
   Array3D<double> Chzey2;

   Array3D<double> Ex_app;
   Array3D<double> Ey_app;
   Array3D<double> Ez_app;

   Array3D<double> Bx_app;
   Array3D<double> By_app;
   Array3D<double> Bz_app;

   Array3D<double> Ex_total;
   Array3D<double> Ey_total;
   Array3D<double> Ez_total;

   Array3D<double> Bx_total;
   Array3D<double> By_total;
   Array3D<double> Bz_total;

   std::vector<GaussianBeam> beams{};
   std::vector<CurrentSource> srcs{};

   Array3D<null_t> empty{}; // for the shits and possibly also some giggles...

   std::unordered_map<std::string, Array3D<double>&> em_map = {
      {"Ex", Ex},
      {"Ey", Ey},
      {"Ez", Ez},
      {"Hx", Hx},
      {"Hy", Hy},
      {"Hz", Hz},
      {"Jx", Jx},
      {"Jy", Jy},
      {"Jz", Jz}
   };
};

using emdata_t = EMData;
} // end namespace tf::electromagnetics


#endif //EM_DATA_HPP
