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
