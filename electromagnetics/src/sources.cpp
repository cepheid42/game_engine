#include "sources.hpp"

#include "program_params.hpp"
#include "constants.hpp"

#include <cmath>
#include <print>

namespace tf::electromagnetics {
  [[nodiscard]] compute_t RickerSource::eval(const compute_t t) const {
    constexpr auto Md = 2.0f;
    const auto alpha = math::SQR(constants::pi * freq * (t - Md / freq));
    return (1.0f - 2.0f * alpha) * std::exp(-alpha);
  }

  [[nodiscard]] compute_t SpatialSource::eval(const compute_t t) const {
    auto result = amplitude;
    for (const auto& src : t_srcs) {
      result *= src->eval(t);
    }
    return result;
  } // end SpatialSource::eval

  [[nodiscard]] compute_t BlackmanHarris::eval(const compute_t t) const {
    if (t > duration) { return 1.0f; }

    const auto c1 = std::cos(bn[0] * omega * t);
    const auto c2 = std::cos(bn[1] * omega * t);
    const auto c3 = std::cos(bn[2] * omega * t);
    return an[0] + (an[1] * c1) + (an[2] * c2) + (an[3] * c3);
  }

  [[nodiscard]] compute_t ContinuousSource::eval(const compute_t t) const {
    if (t < start or t > stop) { return 0.0f; }
    // return std::sin(omega * t - phase);
    return ramp.eval(t) * std::sin(omega * t - phase);
    // return std::sin(omega * t - phase);
  }

  void CurrentSource::apply(const compute_t t) const {
    const auto& [x0, x1, y0, y1, z0, z1] = src.offsets;
    const auto val = src.eval(t);

    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        for (size_t k = z0; k < z1; ++k) {
          (*field)(i, j, k) = val; // todo: only hardsources allowed here
        }
      }
    }
  } // end CurrentSource::apply
} // end namespace tf::electromagnetics
