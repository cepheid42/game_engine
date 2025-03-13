#include "bc_data.hpp"

#include "program_params.hpp"
#include "constants.hpp"
#include "math_utils.hpp"

#include <vector>
#include <algorithm>

namespace tf::electromagnetics {
  PMLData::PMLData(const std::size_t nx_, const std::size_t ny_, const std::size_t nz_,
                   const offset_t& offsetsE_, const offset_t& offsetsH_, const bool hi_side)
  : psiE{nx_, ny_, nz_}, psiH{nx_, ny_, nz_}, offsetsE{offsetsE_}, offsetsH{offsetsH_}
  {
    init_coefficients(hi_side);
  }

  void PMLData::init_coefficients(const bool hi_side) {
    constexpr auto hstep = 1.0 / (2.0 * static_cast<double>(PMLDepth));
    const std::vector<double> d = math::linspace(1.0, 0.0, PMLDepth, false);
    const std::vector<double> dh = math::linspace(1.0 - hstep, -hstep, PMLDepth, false);

    constexpr auto sigma_max = (0.8 * (PMLGrade + 1.0)) / (dx * constants::eta0);

    const auto sigma_d = calculate_sigma(d, sigma_max);
    const auto alpha_d = calculate_alpha(d);
    calculate_coeffs(bE, cE, sigma_d, alpha_d);

    const auto sigma_dh = calculate_sigma(d, sigma_max);
    const auto alpha_dh = calculate_alpha(d);
    calculate_coeffs(bH, cH, sigma_dh, alpha_dh);

    if (hi_side) {
     std::ranges::reverse(bE);
     std::ranges::reverse(cE);
     std::ranges::reverse(bH);
     std::ranges::reverse(cH);
    }
  }

  std::vector<double> PMLData::calculate_sigma(const std::vector<double>& d, const double sigma_max) {
    auto sigma_bc(d);
    for (auto& x: sigma_bc) {
      x = sigma_max * std::pow(x, PMLGrade);
    }
    return sigma_bc;
  }

  std::vector<double> PMLData::calculate_alpha(const std::vector<double>& d) {
    auto alpha_bc(d);
    for (auto& x: alpha_bc) {
      x = PMLAlphaMax * std::pow(1.0 - x, 1.0);
    }
    return alpha_bc;
  }

  void PMLData::calculate_coeffs(coeffs_t& b, coeffs_t& c, const std::vector<double>& sigma, const std::vector<double>& alpha) {
    constexpr auto coef1 = -dt / constants::eps0;

    for (auto i = 0zu; i < PMLDepth; i++) {
      constexpr auto kappa_bc = 1.0;
      b[i] = std::exp(coef1 * ((sigma[i] / kappa_bc) + alpha[i]));
      c[i] = (sigma[i] * (b[i] - 1.0)) / (kappa_bc * (sigma[i] + (kappa_bc * alpha[i])));
    }
  }


  BCData::BCData(const std::size_t nx_, const std::size_t ny_, const std::size_t nz_)
  : x0(PMLDepth, ny_, nz_, {1, PMLDepth, 0, ny_, 0, nz_}, {0, PMLDepth - 1, 0, ny_, 0, nz_}),
    y0(nx_, PMLDepth, nz_, {0, nx_, 1, PMLDepth, 0, nz_}, {0, nx_, 0, PMLDepth - 1, 0, nz_}),
    z0(nx_, ny_, PMLDepth, {0, nx_, 0, ny_, 1, PMLDepth}, {0, nx_, 0, ny_, 0, PMLDepth - 1}),
    x1(PMLDepth, ny_, nz_,
       {nx_ - 1 - PMLDepth, nx_ - 1, 0, ny_, 0, nz_},
       {nx_ - 1 - PMLDepth, nx_ - 1, 0, ny_, 0, nz_}),
    y1(nx_, PMLDepth, nz_,
       {0, nx_, ny_ - 1 - PMLDepth, ny_ - 1, 0, nz_},
       {0, nx_, ny_ - 1 - PMLDepth, ny_ - 1, 0, nz_}),
    z1(nx_, ny_, PMLDepth,
       {0, nx_, 0, ny_, nz_ - 1 - PMLDepth, nz_ - 1},
       {0, nx_, 0, ny_, nz_ - 1 - PMLDepth, nz_ - 1})
  {}
}
