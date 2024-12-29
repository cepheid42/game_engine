//
// Created by cepheid on 12/13/24.
//

#ifndef EM_SOURCES_H
#define EM_SOURCES_H

#include <cmath>

#include "electromagnetics.param"

namespace tf::electromagnetics::sources
{
  template<typename T>
  struct TemporalSource {
    virtual ~TemporalSource() = default;
    [[nodiscard]] virtual T eval(T) const = 0;
  };

  template<typename T>
  struct GaussianSource final : TemporalSource<T> {
    [[nodiscard]] auto eval(const double t) const override {
      constexpr auto tol = 1e-12;
      const auto val = std::exp(-1.0 * std::pow((t - delay) / width, power));
      return val > tol ? val : 0.0;
    }

    double width;
    double power;
    double cutoff;
    double delay;
  };

  template<typename T>
  struct RickerSource final : TemporalSource<T> {
    [[nodiscard]] auto eval(const double q) const override {
      constexpr auto Np = 20.0;
      constexpr auto Md = 2.0;

      const auto alpha = (M_PI * (cfl * q / Np - Md)) * (M_PI * (cfl * q / Np - Md));
      return (1.0 - 2.0 * alpha) * std::exp(-alpha);
    }
  };

  template<typename T>
  struct ContinuousSource final : TemporalSource<T> {
    [[nodiscard]] auto eval(const double t) const override {
      if (t < start_time or t > stop_time) { return 0.0; }
      const auto src_val = 1.0 * std::sin(omega * t - phase);
      return src_val;
    }

    double freq;
    double omega;
    double start_time;
    double stop_time;
    double phase;
    double delay;
  };
  //
} // end namespace tf::electromagnetics::sources
#endif //EM_SOURCES_H
