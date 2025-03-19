#include "em_solver.hpp"
#include "update_functors.hpp"
#include "bc_functors.hpp"

namespace tf::electromagnetics {
  EMSolver::EMSolver(const std::size_t nx, const std::size_t ny, const std::size_t nz, const double cfl, const double dt)
  : emdata(nx, ny, nz, cfl, dt),
    bcdata(this->emdata)
  {}

  void EMSolver::updateE() {
    ex_update(emdata.Ex, emdata.Hz, emdata.Hy, emdata.Jx, emdata.Cexe, emdata.Cexhz, emdata.Cexhy, emdata.Cjx, {0, 0, 1, 1, 1, 1});
    ey_update(emdata.Ey, emdata.Hx, emdata.Hz, emdata.Jy, emdata.Ceye, emdata.Ceyhx, emdata.Ceyhz, emdata.Cjy, {1, 1, 0, 0, 1, 1});
    ez_update(emdata.Ez, emdata.Hy, emdata.Hx, emdata.Jz, emdata.Ceze, emdata.Cezhy, emdata.Cezhx, emdata.Cjz, {1, 1, 1, 1, 0, 0});
  }

  void EMSolver::updateH() {
    hx_update(emdata.Hx, emdata.Ey, emdata.Ez, emdata.empty, emdata.Chxh, emdata.Chxey, emdata.Chxez, emdata.empty, {0, 0, 0, 0, 0, 0});
    hy_update(emdata.Hy, emdata.Ez, emdata.Ex, emdata.empty, emdata.Chyh, emdata.Chyez, emdata.Chyex, emdata.empty, {0, 0, 0, 0, 0, 0});
    hz_update(emdata.Hz, emdata.Ex, emdata.Ey, emdata.empty, emdata.Chzh, emdata.Chzex, emdata.Chzey, emdata.empty, {0, 0, 0, 0, 0, 0});
  }

  // void EMSolver::updateEbcs() {
  //   // Ex_x0();
  //   // Ex_x1();
  //   Ey_x0(emdata.Ey, emdata.Hz, emdata.Ceyhz, emdata.bcdata.);
  //   Ey_x1();
  //   Ez_x0();
  //   Ez_x1();
  //
  //   Ex_y0();
  //   Ex_y1();
  //   // Ey_y0();
  //   // Ey_y1();
  //   Ez_y0();
  //   Ez_y1();
  //
  //   Ex_z0();
  //   Ex_z1();
  //   Ey_z0();
  //   Ey_z1();
  //   // Ez_z0();
  //   // Ez_z1();
  // }
  //
  // void EMSolver::updateHbcs() {
  //   Hx_x0();
  //   Hx_x1();
  //   Hy_x0();
  //   Hy_x1();
  //   Hz_x0();
  //   Hz_x1();
  //
  //   Hx_y0();
  //   Hx_y1();
  //   Hy_y0();
  //   Hy_y1();
  //   Hz_y0();
  //   Hz_y1();
  //
  //   Hx_z0();
  //   Hx_z1();
  //   Hy_z0();
  //   Hy_z1();
  //   Hz_z0();
  //   Hz_z1();
  // }


  void EMSolver::apply_srcs(const double t) const {
    for (const auto& src: emdata.srcs) { // todo: may need explicit loop if I want threads here
      src.apply(t);
    }
  }

  void EMSolver::advance(const double t) {
    updateH();
    // updateHBCs();

    apply_srcs(t);

    updateE();
    // updateEBCS();
  }
} // end namespace tf::electromagnetics