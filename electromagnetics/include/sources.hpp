#ifndef EM_SOURCES_HPP
#define EM_SOURCES_HPP

#include "program_params.hpp"
#include "constants.hpp"
#include "array.hpp"
#include "math_utils.hpp"

#include <memory>
#include <vector>
#include <cassert>

namespace tf::electromagnetics {
  struct TemporalSource {
    virtual ~TemporalSource() = default;
    [[nodiscard]] virtual compute_t eval(compute_t) const = 0;
  };


  struct RickerSource final : TemporalSource {
    explicit RickerSource(const compute_t freq_) : freq(freq_) {}

    [[nodiscard]] compute_t eval(compute_t) const override;

    compute_t freq;
  }; // end struct RickerSource

  struct BlackmanHarris final : TemporalSource {
    /*
     * https://en.wikipedia.org/wiki/Window_function#Blackman%E2%80%93Harris_window
     * Blackman-Harris Window is similar to a Gaussian except having better frequency
     * behaviors to avoid more DC components.
     *
     * cutoff_factor is equivalent to N/2 -> 1.55 / (2*freq) -> 1.55 / (omega/pi) -> (1.55*pi)/omega -> 4.867/omega
     * bn is from 2*pi*t/N -> 2*pi*t / (1.55 / freq) -> 2*pi*freq*t / 1.55 -> omega * t / 1.55,
     *    so the array is {1/1.55, 2/1.55, 3/1.55}
     *
     */
    explicit BlackmanHarris(const compute_t dx_) : omega(static_cast<compute_t>(constants::pi2 * constants::c) / (Nl * dx_)), duration{cutoff_factor / omega} {}

    [[nodiscard]] compute_t eval(compute_t) const override;

    compute_t omega;
    compute_t duration;
    static constexpr compute_t Nl = 25.0_fp;
    static constexpr compute_t cutoff_factor{4.867f};
    static constexpr compute_t an[4]{0.3588_fp, -0.4882_fp, 0.1413_fp, -0.0116_fp};
    static constexpr compute_t bn[3]{0.6452_fp, 1.2903_fp, 1.9355_fp};
  }; // end struct BlackmanHarris


  struct GaussianSource final : TemporalSource {
    GaussianSource(const compute_t width_, const compute_t power_, const compute_t delay_)
    : width(width_), power(power_), delay(delay_)
    {}

    [[nodiscard]] compute_t eval(const compute_t t) const override {
      constexpr auto tol = 1e-15_fp;
      const auto val = std::exp(-1.0_fp * std::pow((t - delay) / width, power));
      return val > tol ? val : 0.0_fp;
    }

    compute_t width;
    compute_t power;
    compute_t delay;
  }; // end struct GaussianSource


  struct ContinuousSource final : TemporalSource {
    explicit ContinuousSource(const compute_t omega_, const compute_t phase_, const compute_t start_, const compute_t stop_, const compute_t dx_)
    : omega(omega_), start(start_), stop(stop_), phase(phase_)//, ramp{dx_}
    {}

    [[nodiscard]] compute_t eval(compute_t) const override;

    compute_t omega;
    compute_t start;
    compute_t stop;
    compute_t phase;
    // BlackmanHarris ramp;
  }; // end struct ContinuousSource

  struct SpatialSource {
    using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;
    using offset_t = std::array<std::size_t, 6>;

    SpatialSource(temporal_vec&& ts, const compute_t amp_, const offset_t& o_)
    : amplitude(amp_),
      offsets{o_},
      t_srcs{std::move(ts)}
    {}

    [[nodiscard]] compute_t eval(compute_t) const;

    compute_t amplitude;
    offset_t offsets;
    temporal_vec t_srcs;
  }; // end struct SpatialSource

  struct CurrentSource {
    using array_t = Array3D<compute_t>;
    CurrentSource(array_t* const f, SpatialSource&& s)
    : field(f),
      src(std::move(s))
    {}

    // void apply(compute_t) const;

    array_t* const field;
    SpatialSource src;
  }; // end struct CurrentSource

  struct GaussianBeam : CurrentSource {
    using array_t = Array3D<compute_t>;
    using offset_t = std::array<std::size_t, 6>;

    GaussianBeam(array_t*,
      compute_t,
      compute_t,
      const vec3<compute_t>&,
      SpatialSource&&);

    void apply(compute_t) const;

    compute_t waist_size;
    vec3<compute_t> waist_pos;
    std::vector<compute_t> coeffs;
  };


} // end namespace tf::electromagnetics

#endif //EM_SOURCES_HPP
