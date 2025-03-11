#ifndef EM_DATA_HPP
#define EM_DATA_HPP

#include <iostream>

#include "em_constants.hpp"
#include "em_array.hpp"
#include "em_sources.hpp"

namespace tf::electromagnetics {
  template<typename ex_t, typename ey_t, typename ez_t, typename hx_t, typename hy_t, typename hz_t>
  struct EMData {
    EMData() = delete;

    explicit EMData(const std::size_t nx, const std::size_t ny, const std::size_t nz, const double cfl, const double dt)
    : Ex(nx - 1, ny, nz), Jx(nx - 1, ny, nz), Cexe(nx - 1, ny, nz), Cexhy(nx - 1, ny, nz), Cexhz(nx - 1, ny, nz), Cjx(nx - 1, ny, nz),
      Ey(nx, ny - 1, nz), Jy(nx, ny - 1, nz), Ceye(nx, ny - 1, nz), Ceyhx(nx, ny - 1, nz), Ceyhz(nx, ny - 1, nz), Cjy(nx, ny - 1, nz),
      Ez(nx, ny, nz - 1), Jz(nx, ny, nz - 1), Ceze(nx, ny, nz - 1), Cezhx(nx, ny, nz - 1), Cezhy(nx, ny, nz - 1), Cjz(nx, ny, nz - 1),
      Hx(nx, ny - 1, nz - 1), Chxh(nx, ny - 1, nz - 1), Chxey(nx, ny - 1, nz - 1), Chxez(nx, ny - 1, nz - 1),
      Hy(nx - 1, ny, nz - 1), Chyh(nx - 1, ny, nz - 1), Chyex(nx - 1, ny, nz - 1), Chyez(nx - 1, ny, nz - 1),
      Hz(nx - 1, ny - 1, nz), Chzh(nx - 1, ny - 1, nz), Chzex(nx - 1, ny - 1, nz), Chzey(nx - 1, ny - 1, nz)
    {
      init_coefficients(cfl, dt);
    }

    void init_coefficients(double, double);

    Array3D<ex_t> Ex;
    Array3D<ex_t> Jx;
    Array3D<ex_t> Cexe;
    Array3D<ex_t> Cexhy;
    Array3D<ex_t> Cexhz;
    Array3D<ex_t> Cjx;

    Array3D<ey_t> Ey;
    Array3D<ey_t> Jy;
    Array3D<ey_t> Ceye;
    Array3D<ey_t> Ceyhx;
    Array3D<ey_t> Ceyhz;
    Array3D<ey_t> Cjy;

    Array3D<ez_t> Ez;
    Array3D<ez_t> Jz;
    Array3D<ez_t> Ceze;
    Array3D<ez_t> Cezhx;
    Array3D<ez_t> Cezhy;
    Array3D<ez_t> Cjz;

    Array3D<hx_t> Hx;
    Array3D<hx_t> Chxh;
    Array3D<hx_t> Chxey;
    Array3D<hx_t> Chxez;

    Array3D<hy_t> Hy;
    Array3D<hy_t> Chyh;
    Array3D<hy_t> Chyex;
    Array3D<hy_t> Chyez;

    Array3D<hz_t> Hz;
    Array3D<hz_t> Chzh;
    Array3D<hz_t> Chzex;
    Array3D<hz_t> Chzey;

    std::vector<CurrentSource> srcs{};

    Array3D<void> empty{}; // for the shits and possibly also some giggles...
  };


  template<typename ex_t, typename ey_t, typename ez_t, typename hx_t, typename hy_t, typename hz_t>
  void EMData<ex_t, ey_t, ez_t, hx_t, hy_t, hz_t>::init_coefficients(const double cfl, const double dt) {
    const auto e_coeff = cfl * constants::eta0;
    const auto h_coeff = cfl / constants::eta0;
    const auto j_coeff = dt / constants::eps0;

    // todo: add loss terms
    Cexe.fill(1.0);
    Ceye.fill(1.0);
    Ceze.fill(1.0);
    Chxh.fill(1.0);
    Chyh.fill(1.0);
    Chzh.fill(1.0);

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
  }
}



#endif //EM_DATA_HPP
