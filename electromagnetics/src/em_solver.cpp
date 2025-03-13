#include "em_solver.hpp"
#include "update_functors.hpp"
#include "bc_functors.hpp"

namespace tf::electromagnetics {
  EMSolver::EMSolver(const std::size_t nx, const std::size_t ny, const std::size_t nz, const double cfl, const double dt)
  : emdata(nx, ny, nz, cfl, dt)
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