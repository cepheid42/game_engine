#include "em_solver.hpp"

#include "program_params.hpp"
#include "update_functors.hpp"
#include "bc_functors.hpp"

namespace tf::electromagnetics {
  EMSolver::EMSolver(const std::size_t nx, const std::size_t ny, const std::size_t nz, const compute_t cfl, const compute_t dt)
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

  void EMSolver::updateEBCs() {
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

  void EMSolver::updateHBCs() {
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

  void EMSolver::apply_srcs(const compute_t t) const {
    for (const auto& src: emdata.srcs) { // todo: may need explicit loop if I want threads here
      src.apply(t);
    }
  }

  void EMSolver::advance(const compute_t t) {
    updateH();
    updateHBCs();

    apply_srcs(t);

    updateE();
    updateEBCs();
  }
} // end namespace tf::electromagnetics