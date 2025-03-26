#ifndef CURRENT_DEPOSITION_HPP
#define CURRENT_DEPOSITION_HPP

#include "math_utils.hpp"

#include "particles.hpp"
#include "em_data.hpp"

#include <array>
#include <cmath>

namespace tf::particles {
  inline vec3<compute_t> normalizedLocation(const vec3<compute_t>& loc) {
    const auto sx = loc[0] / dx;
    const auto sy = loc[1] / dy;
    const auto sz = loc[2] / dz;
    return {sx - std::floor(sx), sy - std::floor(sy), sz - std::floor(sz)};
  }

  inline std::array<compute_t, 3> quad_shapes(compute_t v) {
    const auto shape01 = 0.5f * math::SQR(0.5f - v);
    const auto shape2  = 0.75f - math::SQR(v);
    return {shape01, shape01, shape2};
  }


  struct CurrentDeposition {
    using emdata_t = electromagnetics::EMData;

    static void update(const std::vector<Particle>& particles, const compute_t charge, const emdata_t& emdata) {
      static constexpr auto third = 1.0f / 3.0f;
      static constexpr auto sixth = 1.0f / 6.0f;

      for (std::size_t pid = 0; pid < particles.size(); pid++) {
        for (const auto& segment: split_trajectory()) {
          const auto& [i, j, k] = segment.cid;

          const auto xs0 = quad_shapes(normalizedLocation(segment.p0[0]));
          const auto xs1 = quad_shapes(normalizedLocation(segment.p1[0]));

          const auto ys0 = quad_shapes(normalizedLocation(segment.p0[1]));
          const auto ys1 = quad_shapes(normalizedLocation(segment.p1[1]));

          const auto zs0 = quad_shapes(normalizedLocation(segment.p0[2]));
          const auto zs1 = quad_shapes(normalizedLocation(segment.p1[2]));

          const auto cx = segment.weight * charge / (Ayz * dt);
          const auto cy = segment.weight * charge / (Axz * dt);
          const auto cz = segment.weight * charge / (Axy * dt);

          auto wT = 0.0f;

          auto wm = -(xs1[0] - xs0[0]);
          auto wp = xs1[2] - xs0[2];
          for (std::size_t m = 0; m < 3; m++) {
            for (std::size_t n = 0; n < 3; n++) { // todo: n is fastest changing index, but moves y-index
              wT = cx * (third * (ys0[n] * zs0[m] + ys1[n] * zs1[m]) + sixth * (ys1[n] * zs0[m] + ys0[n] * zs1[m]));
              // todo: need to be concerned about edges of grid here, where j - 1 < 0, j + 2 > Ny,... etc
              emdata.Jx(i, j + n - 1, k + m - 1) += wm * wT;
              emdata.Jx(i + 1, j + n - 1, k + m - 1) += wp * wT;
            }
          }

          wm = -(ys1[0] - ys0[0]);
          wp = ys1[2] - ys0[2];
          for (std::size_t m = 0; m < 3; m++) {
            for (std::size_t n = 0; n < 3; n++) {
              wT = cy * (third * (xs0[n] * zs0[m] + xs1[n] * zs1[m]) + sixth * (xs1[n] * zs0[m] + xs0[n] * zs1[m]));
              // todo: need to be concerned about edges of grid here, where j - 1 < 0, j + 2 > Ny,... etc
              emdata.Jy(i + n - 1, j, k + m - 1) += wm * wT;
              emdata.Jy(i + n - 1, j + 1, k + m - 1) += wp * wT;
            }
          }

        }
      }
    }
  };
}
#endif //CURRENT_DEPOSITION_HPP
