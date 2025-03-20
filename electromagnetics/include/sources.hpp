#ifndef EM_SOURCES_HPP
#define EM_SOURCES_HPP

#include "array.hpp"
#include "math_utils.hpp"

#include <memory>
#include <vector>

namespace tf::electromagnetics {
  struct TemporalSource {
    virtual ~TemporalSource() = default;
    [[nodiscard]] virtual double eval(double) const = 0;
  };


  struct RickerSource final : TemporalSource {
    explicit RickerSource(const double freq_) : freq(freq_) {}

    [[nodiscard]] double eval(double) const override;

    double freq;
  }; // end struct RickerSource

  struct SpatialSource {
    using temporal_vec = std::vector<std::unique_ptr<TemporalSource>>;
    using offset_t = std::array<std::size_t, 6>;

    SpatialSource(temporal_vec&& ts, const double amp_, const offset_t& o_)
    : amplitude(amp_),
      offsets{o_},
      t_srcs{std::move(ts)}
    {}

    [[nodiscard]] double eval(double) const;

    double amplitude;
    offset_t offsets;
    temporal_vec t_srcs;
  }; // end struct SpatialSource

  struct CurrentSource {
    using array_t = Array3D<double>;
    CurrentSource(array_t* const f, SpatialSource&& s)
    : field(f),
      src(std::move(s))
    {}

    void apply(double) const;

    array_t* const field;
    SpatialSource src;
  }; // end struct CurrentSource
} // end namespace tf::electromagnetics

#endif //EM_SOURCES_HPP
