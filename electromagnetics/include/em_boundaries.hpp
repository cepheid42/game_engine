#ifndef EM_BOUNDARIES_HPP
#define EM_BOUNDARIES_HPP

#include <array>
#include <vector>
#include <print>

#include "program_params.hpp"
#include "em_params.hpp"
#include "array.hpp"
#include "constants.hpp"
#include "math_utils.hpp"

namespace tf::electromagnetics {
//  struct PeriodicData {
//    using offset_t = std::array<std::size_t, 6>;
//
//    explicit PeriodicData(std::size_t nInterior_, std::size_t hi_idx_, offset_t offsets_)
//    : nInterior(nInterior_), hi_idx(hi_idx_), offsets(offsets_)
//    {}
//
//    std::size_t nInterior;
//    std::size_t hi_idx;
//    offset_t offsets;
//  };

  struct PMLData {
    using offset_t = std::array<std::size_t, 6>;
    using coeffs_t = std::array<double, PMLDepth>;

    explicit PMLData(std::size_t nx_, std::size_t ny_, std::size_t nz_, offset_t offsetsE_, offset_t offsetsH_)
    : psiE(nx_, ny_, nz_), psiH(nx_, ny_, nz_), offsetsE(offsetsE_), offsetsH(offsetsH_)
    {
      init_coefficients();
    }

    void init_coefficients() {
      constexpr auto coef1 = -dt / constants::eps0;

      auto hstep = 1.0 / (2.0 * static_cast<double>(PMLDepth));
      std::vector<double> d = math::linspace(1.0, 0.0, PMLDepth, false);
      std::vector<double> dh = math::linspace(1.0 - hstep, -hstep, PMLDepth, false);

	  constexpr auto sigma_max = (0.8 * (PMLGrade + 1.0)) / (dx * constants::eta0);

	  auto sigma_bc(d);
	  auto alpha_bc(d);
	  auto kappa_bc(d.size(), 1.0);

	  for (auto& x: sigma_bc) {
		x = sigma_max * std::pow(x, PMLGrade);
	  }

	  for (auto& x: alpha_bc) {
		x = PMLAlphaMax * std::pow(1.0 - x, 1.0);
      }

	  for (auto i = 0zu; i < PMLDepth; i++) {
		bE[i] = std::exp(coef1 * ((sigma_bc[i] / kappa_bc[i]) + alpha_bc[i]));
		cE[i] = (sigma_bc[i] * (b[i] - 1.0)) / (kappa_bc[i] * (sigma_bc[i] + (kappa_bc[i] * alpha_bc[i])));
	  }

		// todo: add H coeffs
		// todo: how to make reverse option for hi side?


    }

    Array3D<double> psiE;
    Array3D<double> psiH;
    offset_t offsetsE;
    offset_t offsetsH;
    coeffs_t bE{};
    coeffs_t cE{};
    coeffs_t bH{};
    coeffs_t cH{};
  };


  struct BCData {
    BCData(std::size_t nx_, std::size_t ny_, std::size_t nz_)
    : x0(PMLDepth, ny_, nz_, {1, PMLDepth, 0,      ny_, 0,      nz_}, {0, PMLDepth - 1, 0,          ny_, 0,          nz_}),
      y0(nx_, PMLDepth, nz_, {0,      nx_, 1, PMLDepth, 0,      nz_}, {0,          nx_, 0, PMLDepth - 1, 0,          nz_}),
      z0(nx_, ny_, PMLDepth, {0,      nx_, 0,      ny_, 1, PMLDepth}, {0,          nx_, 0,          ny_, 0, PMLDepth - 1}),
      x1(PMLDepth, ny_, nz_, {nx_ - 1 - PMLDepth, nx_ - 1,                  0,     ny_,                   0,     nz_}, {nx_ - 1 - PMLDepth, nx_ - 1,                  0,     ny_,                  0,     nz_}),
      y1(nx_, PMLDepth, nz_, {                 0,     nx_, ny_ - 1 - PMLDepth, ny_ - 1,                   0,     nz_}, {                 0,     nx_, ny_ - 1 - PMLDepth, ny_ - 1,                  0,     nz_}),
      z1(nx_, ny_, PMLDepth, {                 0,     nx_,                  0,      ny_, nz_ - 1 - PMLDepth, nz_ - 1}, {                 0,     nx_,                  0,     ny_, nz_ - 1 - PMLDepth, nz_ - 1})
    {}

    PMLData x0;
    PMLData y0;
    PMLData z0;
    PMLData x1;
    PMLData y1;
    PMLData z1;
  };
}

#endif //EM_BOUNDARIES_HPP
