#ifndef EM_SOURCES_HPP
#define EM_SOURCES_HPP

#include "program_params.hpp"
#include "constants.hpp"
#include "array.hpp"
#include "math_utils.hpp"

#include <memory>
#include <vector>

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
    explicit BlackmanHarris(const compute_t dx_) : omega(constants::pi2 * constants::c / (Nl * dx_)), duration{cutoff_factor / omega} {}

    [[nodiscard]] compute_t eval(compute_t) const override;

    compute_t omega;
    compute_t duration;
    static constexpr compute_t Nl = 25.0f;
    static constexpr compute_t cutoff_factor{4.867f};
    static constexpr compute_t an[4]{0.3588f, -0.4882f, 0.1413f, -0.0116f};
    static constexpr compute_t bn[3]{0.6452f, 1.2903f, 1.9355f};
  }; // end struct BlackmanHarris
  
  struct ContinuousSource final : TemporalSource {
    explicit ContinuousSource(const compute_t omega_, const compute_t phase_, const compute_t start_, const compute_t stop_, const compute_t dx_)
    : omega(omega_), start(start_), stop(stop_), phase(phase_), ramp{dx_}
    {}

    [[nodiscard]] compute_t eval(compute_t) const override;

    compute_t omega;
    compute_t start;
    compute_t stop;
    compute_t phase;
    BlackmanHarris ramp;
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

    void apply(compute_t) const;

    array_t* const field;
    SpatialSource src;
  }; // end struct CurrentSource
} // end namespace tf::electromagnetics

#endif //EM_SOURCES_HPP
