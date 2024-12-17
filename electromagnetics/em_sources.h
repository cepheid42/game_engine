//
// Created by cepheid on 12/13/24.
//

#ifndef EM_SOURCES_H
#define EM_SOURCES_H

namespace tf::electromagnetics::sources
{
  template<typename Derived>
  struct TemporalSource {
    [[nodiscard]] auto eval(const double t) const {
      return static_cast<Derived*>(this)->eval(t);
    }
  };

  struct GaussianSource : TemporalSource<GaussianSource> {
    [[nodiscard]] auto eval(const double t) const {
      constexpr auto tol = 1e-12;
      const auto val = std::exp(-1.0 * std::pow((t - delay) / width, power));
      return val > tol ? val : 0.0;
    }

    double width;
    double power;
    double cutoff;
    double delay;
  };

  struct ContinuousSource : TemporalSource<ContinuousSource> {
    [[nodiscard]] auto eval(const double t) const {
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
} // end namespace tf::electromagnetics::sources
#endif //EM_SOURCES_H
