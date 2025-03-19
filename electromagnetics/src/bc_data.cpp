#include "bc_data.hpp"

#include "program_params.hpp"
#include "constants.hpp"
#include "math_utils.hpp"

#include <vector>
#include <algorithm>

namespace tf::electromagnetics {
  // PMLData::PMLData(const std::size_t nx_, const std::size_t ny_, const std::size_t nz_,
  //                  const offset_t& offsets_, const bool isHi, const bool isH)
  // : psi{nx_, ny_, nz_}, offsets{offsets_}
  // {
  //   init_coefficients(isHi, isH);
  // }
  //
  // void PMLData::init_coefficients(const bool isHi, const bool isHField) {
  //   std::vector<double> d = math::linspace(1.0, 0.0, PMLDepth, false);
  //   // const std::vector<double> dh = math::linspace(1.0 - hstep, -hstep, PMLDepth, false);
  //
  //   if constexpr (isHField) {
  //     constexpr auto hstep = 1.0 / (2.0 * static_cast<double>(PMLDepth));
  //     for (auto& x: d) { x -= hstep; }
  //   }
  //
  //   constexpr auto sigma_max = (0.8 * (PMLGrade + 1.0)) / (dx * constants::eta0);
  //
  //   const auto sigma_d = calculate_sigma(d, sigma_max);
  //   const auto alpha_d = calculate_alpha(d);
  //   calculate_coeffs(b, c, sigma_d, alpha_d);
  //
  //   if constexpr (isHi) {
  //    std::ranges::reverse(b);
  //    std::ranges::reverse(c);
  //   }
  // }
  //
  // std::vector<double> PMLData::calculate_sigma(const std::vector<double>& d, const double sigma_max) {
  //   auto sigma_bc(d);
  //   for (auto& x: sigma_bc) {
  //     x = sigma_max * std::pow(x, PMLGrade);
  //   }
  //   return sigma_bc;
  // }
  //
  // std::vector<double> PMLData::calculate_alpha(const std::vector<double>& d) {
  //   auto alpha_bc(d);
  //   for (auto& x: alpha_bc) {
  //     x = PMLAlphaMax * std::pow(1.0 - x, 1.0);
  //   }
  //   return alpha_bc;
  // }
  //
  // void PMLData::calculate_coeffs(coeffs_t& b, coeffs_t& c, const std::vector<double>& sigma, const std::vector<double>& alpha) {
  //   constexpr auto coef1 = -dt / constants::eps0;
  //
  //   for (auto i = 0zu; i < PMLDepth; i++) {
  //     constexpr auto kappa_bc = 1.0;
  //     b[i] = std::exp(coef1 * ((sigma[i] / kappa_bc) + alpha[i]));
  //     c[i] = (sigma[i] * (b[i] - 1.0)) / (kappa_bc * (sigma[i] + (kappa_bc * alpha[i])));
  //   }
  // }
}
