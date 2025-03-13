#include "constants.hpp"
#include "sources.hpp"

#include <cmath>

namespace tf::electromagnetics {
  [[nodiscard]] double RickerSource::eval(const double t) const {
    constexpr auto Md = 2.0;
    const auto alpha = SQR(constants::pi * freq * (t - Md / freq));
    return (1.0 - 2.0 * alpha) * std::exp(-alpha);
  }

  [[nodiscard]] double SpatialSource::eval(const double t) const {
    auto result = amplitude;
    for (const auto& src : t_srcs) {
      result *= src->eval(t);
    }
    return result;
  } // end SpatialSource::eval

  void CurrentSource::apply(const double t) const {
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
