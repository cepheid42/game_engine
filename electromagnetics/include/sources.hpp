#ifndef EM_SOURCES_HPP
#define EM_SOURCES_HPP

#include "program_params.hpp"
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
