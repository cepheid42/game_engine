#include "em_data.hpp"

#include "program_params.hpp"
#include "constants.hpp"

#include <print>

namespace tf::electromagnetics {
  EMData::EMData(const std::size_t nx, const std::size_t ny, const std::size_t nz, const compute_t cfl, const compute_t dt)
  : Ex(nx - 1, ny, nz), Jx(nx - 1, ny, nz), Cexe(nx - 1, ny, nz), Cexhy(nx - 1, ny, nz), Cexhz(nx - 1, ny, nz), Cjx(nx - 1, ny, nz),
    Ey(nx, ny - 1, nz), Jy(nx, ny - 1, nz), Ceye(nx, ny - 1, nz), Ceyhx(nx, ny - 1, nz), Ceyhz(nx, ny - 1, nz), Cjy(nx, ny - 1, nz),
    Ez(nx, ny, nz - 1), Jz(nx, ny, nz - 1), Ceze(nx, ny, nz - 1), Cezhx(nx, ny, nz - 1), Cezhy(nx, ny, nz - 1), Cjz(nx, ny, nz - 1),
    Hx(nx, ny - 1, nz - 1), Chxh(nx, ny - 1, nz - 1), Chxey(nx, ny - 1, nz - 1), Chxez(nx, ny - 1, nz - 1),
    Hy(nx - 1, ny, nz - 1), Chyh(nx - 1, ny, nz - 1), Chyex(nx - 1, ny, nz - 1), Chyez(nx - 1, ny, nz - 1),
    Hz(nx - 1, ny - 1, nz), Chzh(nx - 1, ny - 1, nz), Chzex(nx - 1, ny - 1, nz), Chzey(nx - 1, ny - 1, nz)
  {
    init_coefficients(cfl, dt);
  }

  void EMData::init_coefficients(const compute_t cfl, const compute_t dt) {
    const auto e_coeff = cfl * constants::eta0;
    const auto h_coeff = cfl / constants::eta0;
    const auto j_coeff = dt / constants::eps0;

    // todo: add loss terms
    Cexe.fill(1.0f);
    Ceye.fill(1.0f);
    Ceze.fill(1.0f);
    Chxh.fill(1.0f);
    Chyh.fill(1.0f);
    Chzh.fill(1.0f);

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

  compute_t EMData::Bx(const std::size_t i, const std::size_t j, const std::size_t k) const { return constants::mu0 * Hx(i, j, k); }
  compute_t EMData::By(const std::size_t i, const std::size_t j, const std::size_t k) const { return constants::mu0 * Hy(i, j, k); }
  compute_t EMData::Bz(const std::size_t i, const std::size_t j, const std::size_t k) const { return constants::mu0 * Hz(i, j, k); }



} // end namespace tf::electromagnetics

