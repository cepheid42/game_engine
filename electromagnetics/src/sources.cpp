#include "sources.hpp"

#include "program_params.hpp"
#include "constants.hpp"

#include <cmath>
#include <print>

namespace tf::electromagnetics {
  [[nodiscard]] compute_t RickerSource::eval(const compute_t t) const {
    constexpr auto Md = 2.0_fp;
    const auto alpha = math::SQR(static_cast<compute_t>(constants::pi<compute_t>) * freq * (t - Md / freq));
    return (1.0_fp - 2.0_fp * alpha) * std::exp(-alpha);
  }

  [[nodiscard]] compute_t BlackmanHarris::eval(const compute_t t) const {
    if (t > duration) { return 1.0_fp; }

    const auto c1 = std::cos(bn[0] * omega * t);
    const auto c2 = std::cos(bn[1] * omega * t);
    const auto c3 = std::cos(bn[2] * omega * t);
    return an[0] + (an[1] * c1) + (an[2] * c2) + (an[3] * c3);
  }

  [[nodiscard]] compute_t GaussianSource::eval(const compute_t t) const {
    constexpr auto tol = 1e-15_fp;
    const auto val = std::exp(-1.0_fp * std::pow((t - delay) / width, power));
    // std::println("{}", val);
    return val > tol ? val : 0.0_fp;
  }

  [[nodiscard]] compute_t ContinuousSource::eval(const compute_t t) const {
    if (t < start or t > stop) { return 0.0_fp; }
    const auto val = std::sin(omega * t - phase);
    // std::println("{}", val);
    return val;

    // return ramp.eval(t) * std::sin(omega * t - phase);
  }

  [[nodiscard]] compute_t SpatialSource::eval(const compute_t t) const {
    auto result = amplitude;
    for (const auto& src : t_srcs) {
      result *= src->eval(t);
    }
    return result;
  } // end SpatialSource::eval

  // void CurrentSource::apply(const compute_t t) const {
  //   const auto& [x0, x1, y0, y1, z0, z1] = src.offsets;
  //   const auto val = src.eval(t);
  //
  //   for (size_t i = x0; i < x1; ++i) {
  //     for (size_t j = y0; j < y1; ++j) {
  //       for (size_t k = z0; k < z1; ++k) {
  //         (*field)(i, j, k) = val; // todo: only hardsources allowed here
  //       }
  //     }
  //   }
  // } // end CurrentSource::apply

  GaussianBeam::GaussianBeam(array_t* const f_,
                             const compute_t waist_,
                             const compute_t omega_,
                             const vec3<compute_t>& waist_pos_,
                             SpatialSource&& s_)
  : CurrentSource(f_, std::forward<SpatialSource>(s_)),
    waist_size(waist_),
    waist_pos(waist_pos_),
    coeffs(src.offsets[5] - src.offsets[4])
  {
    const auto& [x0, x1, y0, y1, z0, z1] = src.offsets;

    const auto xpos = x_range[0] + (static_cast<compute_t>(x0) * dx);
    // const auto z = waist_pos[0] - xpos; // +x direction
    const auto z = xpos - waist_pos[0]; // -x direction

    assert(z != 0.0_fp);
    const auto k = omega_ / static_cast<compute_t>(constants::c<compute_t>);
    const auto zR = 0.5_fp * k * math::SQR(waist_size);
    const auto wz = waist_size * std::sqrt(1.0_fp + math::SQR(z / zR));

    const auto RC = z * (1.0_fp + math::SQR(zR / z));
    const auto gouy = std::atan(z / zR);
    const auto c1 = waist_size / wz;

    // std::println("k: {}\nz: {}\nzR: {}\nw(z): {}\nRC: {}\nGouy: {}\nc1: {}", k, z, zR, wz, RC, gouy, c1);

    std::vector<compute_t> r(z1 - z0);

    for (std::size_t i = z0; i < z1; ++i) {
      r[i - z0] = z_range[0] + dz * static_cast<compute_t>(i);
    }

    for (std::size_t i = 0; i < r.size(); ++i) {
      // std::println("{}, {}",std::exp(-1.0_fp * math::SQR(r[i] / wz)), std::cos((k * z) + (k * math::SQR(r[i]) / (2.0_fp * RC)) - gouy));
      coeffs[i] = c1 * std::exp(-1.0_fp * math::SQR(r[i] / wz)) * std::cos((k * z) + (k * math::SQR(r[i]) / (2.0_fp * RC)) - gouy);
    }
  } // end GaussianBeam ctor

  void GaussianBeam::apply(const compute_t t) const {
    const auto& [x0, x1, y0, y1, z0, z1] = src.offsets;
    const auto val = src.eval(t);
    for (size_t i = x0; i < x1; ++i) {
      for (size_t j = y0; j < y1; ++j) {
        for (size_t k = z0; k < z1; ++k) {
          (*field)(i, j, k) = coeffs[k - z0] * val;
          // (*field)(i, j, k) = val;
        }
      }
    }
  }
} // end namespace tf::electromagnetics
