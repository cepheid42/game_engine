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
    GaussianSource(T width_, T power_, T delay_)
    : width(width_), power(power_), delay(delay_)
    {}

    [[nodiscard]] T eval(const T t) const override {
      // constexpr auto tol = 1e-12;
      const auto val = std::exp(-1.0 * std::pow((t - delay) / width, power));
      // DBG(val);
      return val;
    }

    T width;
    T power;
    // T cutoff;
    T delay;
  };

  template<typename T>
  struct RickerSource final : TemporalSource<T> {
    explicit RickerSource(const T freq_) : freq(freq_) {}

    [[nodiscard]] T eval(const T t) const override {
      constexpr auto Md = 2.0;

      const auto alpha = SQR(M_PI * freq * (t - Md / freq));
      return (1.0 - 2.0 * alpha) * std::exp(-alpha);
    }

    T freq;
  };

  template<typename T>
  struct ContinuousSource final : TemporalSource<T> {
    ContinuousSource(T omega_, T start_, T stop_, T phase_)
    : omega(omega_), start(start_), stop(stop_), phase(phase_)
    {}

    [[nodiscard]] T eval(const T t) const override {
      if (t < start or t > stop) { return 0.0; }
      // DBG(std::sin(omega * t - phase));
      return std::sin(omega * t - phase);
    }

    // T freq;
    T omega;
    T start;
    T stop;
    T phase;
    // T delay;
  };


  template<typename T>
  struct SpatialSource {
    SpatialSource(const T amp_, const size_t x0_, const size_t x1_, const size_t y0_=0, const size_t y1_=1, const size_t z0_=0, const size_t z1_=1)
    : amp(amp_), x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), t_srcs{}
    {}

    [[nodiscard]] T eval(const T t) {
      auto result = amp;
      for (const auto& src : t_srcs) {
        result *= src->eval(t);
      }

      // DBG(result);
      return result;
    }

    T amp;
    size_t x0, x1, y0, y1, z0, z1;
    std::vector<std::unique_ptr<TemporalSource<T>>> t_srcs;
  };

  // todo: Should I fold this into the SpatialSource class to get rid of a level of abstraction?
  template<typename Array>
  struct CurrentSource {
    using array_t = Array;
    using value_t = typename Array::value_t;
    using dimension_t = typename Array::dimension_t;

    CurrentSource(Array* f, SpatialSource<value_t>&& s)
    : field(f),
      src(std::move(s))
    {}

    // todo: These are soft sources, but don't fix the soft source problem (incorrect amplitude)
    //       Not sure what the best fix is. Auxilliary sources seem overkill. Maybe a lookup table?
    void apply(const value_t q) const
    requires (dimension_t::value == 1)
    {
      for (size_t i = src.x0; i <= src.x1; ++i) {
        (*field)[i] += src.eval(q);
      }
    }

    void apply(const value_t q)
    requires (dimension_t::value == 2)
    {
      for (size_t i = src.x0; i < src.x1; ++i) {
        for (size_t j = src.y0; j < src.y1; ++j) {
          (*field)(i, j) += src.eval(q);
        }
      }
    }

    void apply(const value_t q) const
    requires (dimension_t::value == 3)
    {
      for (size_t i = src.x0; i < src.x1; ++i) {
        for (size_t j = src.y0; j < src.y1; ++j) {
          for (size_t k = src.z0; k < src.z1; ++k) {
            (*field)(i, j, k) += src.eval(q);
          }
        }
      }
    }

    Array* const field;
    SpatialSource<value_t> src;
  };

  template<typename Array>
  struct GaussianBeam : SpatialSource<typename Array::value_t> {
    using value_t = typename Array::value_t;
    using vec2_t = std::array<value_t, 2>;

    GaussianBeam(Array* const f,
                 const value_t amp_,
                 const value_t w0_,
                 const value_t omega_,
                 const vec2_t& p0_,
                 const size_t x0_,
                 const size_t x1_,
                 const size_t y0_,
                 const size_t y1_,
                 const value_t dx)
    : SpatialSource<value_t>{amp_, x0_, x1_, y0_, y1_},
      field(f), w0(w0_), omega(omega_), p0{p0_}, Ecoeffs(y1_ - y0_)
    {
      constexpr auto c0 = 299792458.0;

      const auto z = (5.0 * dx) - p0[0]; // +x direction

      const auto k = omega / c0;
      const auto z_R = 0.5 * k * SQR(w0);
      const auto w_z = w0 * sqrt(1.0 + SQR(z / z_R));

      const auto RC = z * (1.0 + SQR(z_R / z));
      const auto gouy = std::atan(z / z_R);
      const auto c1 = (w0 / w_z);

      constexpr auto offset = dx * 5.0; // todo: number of nodes the value_tFSF is inset by, how to generalize this?
      std::vector<double> r(y1_ - y0_, offset);

      for (size_t i = 0; i < r.size(); ++i) {
        r[i] += dx * static_cast<double>(i) - 0.5; // todo: 0.5 is the center of the y range
      }

      for (size_t i = 0; i < r.size(); ++i) {
        Ecoeffs[i] = c1 * std::exp(-1.0 * SQR(r[i] / w_z)) * std::cos((k * z) + ((k * SQR(r[i])) / (2.0 * RC)) - gouy);
      }
    }

    void apply(const value_t q)
    {
      // todo: Soft source? Hard source? Who know! Find out next time on... Sourcing Sources with Sorcery!
      for (size_t i = SpatialSource<value_t>::x0; i < SpatialSource<value_t>::x1; ++i) {
        for (size_t j = SpatialSource<value_t>::y0; j < SpatialSource<value_t>::y1; ++j) {
          const auto Eidx = j - SpatialSource<value_t>::y0;
          (*field)(i, j) += Ecoeffs[Eidx] * SpatialSource<value_t>::eval(q);
        }
      }
    }

    Array* const field;
    value_t w0;
    value_t omega;
    vec2_t p0;
    std::vector<value_t> Ecoeffs;
  };
  
  //
} // end namespace tf::electromagnetics::sources
#endif //EM_SOURCES_H
